"""
04 degree distribution visual.py
================================
Visualizes the degree distribution of networks with M_k annotations.

This script computes and plots the mean degree of "common friends to k" (M_k)
for k = 0, 1, 2, ... overlaid on the degree distribution histogram.

Key concepts from the paper:
- M_0: Mean degree in the population
- M_1 (M_f): Mean degree of friends (the classic friendship paradox)
- M_2: Mean degree of common friends to 2 people
- M_k: Mean degree of common friends to k people

The central finding is that M_k increases with k — common friends to more
people tend to have dramatically higher degree.

Mathematical framework (see Appendix B of the paper):
- F_{d,k} = d(d-1)(d-2)...(d-k+1) is the "falling factorial"
  This counts the number of ways someone with degree d can be a common friend to k people.
- M_k = Σ d_i × F_{d_i,k} / Σ F_{d_i,k}  (weighted mean, Equation 7)
- V_k = Σ (d_i - M_k)² × F_{d_i,k} / Σ F_{d_i,k}  (weighted variance, Equation 8)
"""

from common import *

# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_degree_distribution_with_Mk(
    degrees,
    log_scale=True,
    degree_cutoff=0,
    show_Mk_up_to=5,
    bin_width=1,
    m0_label_lift=0,
    tail_tick_top_n=200,
    tail_tick_threshold=None,
):
    """
    Plot degree distribution histogram with M_k values annotated.
    
    Creates a histogram of the degree distribution and overlays vertical red
    lines at M_0, M_1, M_2, ... showing how common friends to k people have
    progressively higher mean degree.
    
    Args:
        degrees: NumPy array of node degrees
        log_scale: If True, use logarithmic y-axis (useful for heavy-tailed distributions)
        degree_cutoff: Minimum degree to include in histogram
        show_Mk_up_to: Show M_0 through M_{show_Mk_up_to}, or a list/range of specific k values
        bin_width: Width of histogram bins
        m0_label_lift: Extra vertical lift for the M_0 label (to avoid overlap)
        tail_tick_top_n: Number of highest-degree nodes drawn as individual tail
            tick marks. Default 100. Ignored if ``tail_tick_threshold`` is given.
        tail_tick_threshold: Absolute degree above which to draw tail ticks.
            Overrides ``tail_tick_top_n`` when not None.
    
    Returns:
        (M_dict, V_dict): The computed M_k and V_k values
    """
    M, V = compute_Mk_and_Vk(degrees)
    print_mk(M)

    # PNAS single-column width (8.7 cm ≈ 3.42 in); height scaled so aspect is
    # comparable to the old two-column panels without wasting vertical space.
    _fs_mk = 10
    plt.figure(figsize=(PNAS_WIDTHS_IN["1col"], 2.45))
    
    # Plot histogram. Note: range() is exclusive of stop, so we add bin_width
    # to include the maximum-degree node (otherwise small networks lose the tail).
    degrees_to_plot = degrees[degrees >= degree_cutoff]
    max_degree = int(max(degrees))
    plt.hist(
        degrees_to_plot,
        bins=range(degree_cutoff, max_degree + bin_width + 1, bin_width),
        color='gray',
        align='mid'
    )
    
    # Add vertical line at degree=1 (common in social networks)
    count_degree_1 = np.sum(degrees == 1)
    if count_degree_1 > 0:
        plt.vlines(1, 0, count_degree_1, color='gray', alpha=1)
    
    if log_scale:
        plt.yscale('log')
    
    # Adjust axis limits
    x_min, x_max = plt.xlim()
    plt.xlim(x_min, x_max)
    plt.ylim(1, None)
    
    y_min, y_max = plt.ylim()
    y_range = y_max - y_min
    
    # Compute annotation heights (different for log vs linear scale)
    if log_scale:
        height_tick_marks = y_range ** 0.1
        height_Mk_labels = y_range ** 0.4
        height_M0_label = y_range ** (0.4 + m0_label_lift)
    else:
        height_tick_marks = y_range * 0.1
        height_Mk_labels = y_range * 0.4
        height_M0_label = y_range * (0.4 + m0_label_lift)
    
    if show_Mk_up_to:
        # Tail ticks: explicit absolute threshold wins; otherwise draw the
        # `tail_tick_top_n` highest-degree nodes (ties at the cutoff included).
        if tail_tick_threshold is not None:
            high_degree_nodes = degrees[(degrees > tail_tick_threshold) & (degrees <= x_max)]
        elif tail_tick_top_n and tail_tick_top_n > 0:
            in_range = degrees[degrees <= x_max]
            if in_range.size <= tail_tick_top_n:
                high_degree_nodes = in_range
            else:
                # Cutoff = N-th largest degree; include ties so we don't drop
                # arbitrary nodes at the boundary (may yield slightly > N ticks).
                cutoff = float(np.partition(in_range, -tail_tick_top_n)[-tail_tick_top_n])
                high_degree_nodes = in_range[in_range >= cutoff]
        else:
            high_degree_nodes = np.empty(0, dtype=degrees.dtype)
        plt.vlines(high_degree_nodes, 0, height_tick_marks, color='black', alpha=0.2)
        
        # Draw labeled vertical lines for each M_k
        for k, m_k in M.items():
            # Check if we should show this k
            if isinstance(show_Mk_up_to, (list, set, range)):
                if k not in show_Mk_up_to:
                    continue
            elif isinstance(show_Mk_up_to, int):
                if k > show_Mk_up_to:
                    break
            
            # Create label: M for k=0, M_f for k=1, M_k otherwise
            if k == 0:
                label = "$M$"
            elif k == 1:
                label = "$M_f$"
            else:
                label = f"$M_{k}$"
            
            # Draw line and label (M_0 gets extra height to stand out).
            # One step above axis tick size so M_k stays readable without dominating.
            if k == 0:
                plt.vlines(m_k, 0, height_M0_label, color='red', alpha=0.5)
                plt.text(m_k, height_M0_label, label, color='red',
                         ha='center', va='bottom', fontsize=_fs_mk)
            else:
                plt.vlines(m_k, 0, height_Mk_labels, color='red', alpha=0.5)
                plt.text(m_k, height_Mk_labels, label, color='red',
                         ha='center', va='bottom', fontsize=_fs_mk)
    
    # Format y-axis ticks for log scale
    if log_scale:
        max_power = int(np.log10(y_max)) + 1
        y_ticks = [10 ** i for i in range(max_power)]
        plt.yticks(y_ticks, [f'{tick:,.0f}' for tick in y_ticks])

    return M, V


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    _fs_axes = 9
    _fs_ticks = 8
    for net in NETWORKS:
        title = {
            "neworleans": "NEW ORLEANS FACEBOOK NETWORK",
            "citations": "CITATION NETWORK",
            "indian_microfinance": "INDIAN MICROFINANCE VILLAGE NETWORK",
            "power": "WESTERN US POWER GRID",
            "highschool": "MARSEILLE HIGH SCHOOL 2013 PROXIMITY NETWORK",
            "celegans_chemical": "C. ELEGANS CHEMICAL CONNECTOME",
        }.get(net.name, net.name.upper())

        print("\n" + "=" * 60)
        print(title)
        print("=" * 60)

        # Tail ticks default to the 100 highest-degree nodes; works across
        # heavy-tailed and small networks without per-network tuning.
        plot_degree_distribution_with_Mk(
            net.ds,
            log_scale=True,
            show_Mk_up_to=range(0, 6) if net.name == "citations" else 5,
            bin_width=1,
            m0_label_lift=0.1 if net.name == "citations" else 0,
            tail_tick_top_n=100,
        )
        if net.name == "citations":
            plt.xlabel("In-degree (number of citations)", fontsize=_fs_axes)
        else:
            plt.xlabel("Degree", fontsize=_fs_axes)
        plt.ylabel("Count", fontsize=_fs_axes)
        plt.tick_params(axis='both', labelsize=_fs_ticks)
        written = save_figure(f"{net.name}.degree-withM")
        plt.close()
        for path in written:
            print(f"Saved: {path}")
