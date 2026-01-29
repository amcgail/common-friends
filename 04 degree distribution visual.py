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
# FALLING FACTORIAL COMPUTATION
# =============================================================================

_falling_factorial_cache = {}

def falling_factorial(d, k):
    """
    Compute the falling factorial F_{d,k} = d × (d-1) × (d-2) × ... × (d-k+1).
    
    This counts the number of ordered k-tuples of a node's friends, i.e.,
    the number of ways a node with degree d can appear as a common friend
    to k specific people.
    
    Examples:
        falling_factorial(5, 2) = 5 × 4 = 20
        falling_factorial(5, 3) = 5 × 4 × 3 = 60
        falling_factorial(3, 3) = 3 × 2 × 1 = 6
        falling_factorial(2, 3) = 0  (can't be common friend to 3 with only 2 friends)
    
    Args:
        d: Degree of the node
        k: Number of common friends required
    
    Returns:
        The falling factorial F_{d,k}, or 0 if d < k
    """
    if (d, k) in _falling_factorial_cache:
        return _falling_factorial_cache[d, k]
    
    result = 1
    for i in range(k):
        result *= (d - i)
    
    _falling_factorial_cache[d, k] = result
    return result


def falling_factorial_array(max_degree, k):
    """
    Efficiently compute falling factorials for degrees k, k+1, ..., max_degree.
    
    Uses a recurrence relation for efficiency:
        F_{d+1,k} = F_{d,k} × (d+1) / (d+1-k)
    
    This avoids recomputing the full product for each degree.
    
    Args:
        max_degree: Maximum degree to compute
        k: The k in F_{d,k}
    
    Returns:
        NumPy array where array[d-k] = F_{d,k} for d = k, k+1, ..., max_degree
    """
    # Start with F_{k,k} = k!
    F_k_k = falling_factorial(k, k)
    weights = [F_k_k]
    
    # Build up using recurrence: F_{d+1,k} = F_{d,k} × (d+1) / (d+1-k)
    for d in range(k + 1, max_degree + 1):
        next_weight = weights[-1] * d / (d - k)
        weights.append(next_weight)
    
    return np.array(weights)


# =============================================================================
# M_k AND V_k COMPUTATION
# =============================================================================

def compute_Mk_and_Vk(degrees, max_k=25):
    """
    Compute mean degree (M_k) and variance (V_k) for common friends to k people.
    
    For each k from 0 to max_k, this computes:
    - M_k: The weighted mean degree, where each node is weighted by F_{d,k}
    - V_k: The weighted variance of degree
    
    These implement Equations 7-8 from the paper's Appendix B:
        M_k = Σ d_i × F_{d_i,k} / Σ F_{d_i,k}
        V_k = Σ (d_i - M_k)² × F_{d_i,k} / Σ F_{d_i,k}
    
    Args:
        degrees: NumPy array of node degrees
        max_k: Maximum k to compute (default 25)
    
    Returns:
        (M_dict, V_dict): Dictionaries mapping k -> M_k and k -> V_k
    """
    M = {0: degrees.mean()}
    V = {0: degrees.var()}
    
    for k in range(1, max_k):
        # Compute falling factorial weights for all possible degrees
        weights_by_degree = falling_factorial_array(max(degrees), k)
        
        # Filter to nodes with degree >= k (others have F_{d,k} = 0)
        eligible_degrees = degrees[degrees >= k]
        
        if len(eligible_degrees) == 0:
            M[k] = np.nan
            V[k] = np.nan
            continue
        
        # Get weight for each eligible node: weights_by_degree[d - k] = F_{d,k}
        node_weights = weights_by_degree[eligible_degrees - k]
        
        # Weighted mean: M_k = Σ d × w / Σ w
        M[k] = (node_weights * eligible_degrees).sum() / node_weights.sum()
        
        # Weighted variance: V_k = Σ (d - M_k)² × w / Σ w
        V[k] = (node_weights * (eligible_degrees - M[k])**2).sum() / node_weights.sum()
    
    return M, V


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_degree_distribution_with_Mk(
    degrees,
    log_scale=True,
    degree_cutoff=0,
    show_Mk_up_to=5,
    bin_width=1,
    m0_label_lift=0
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
    
    Returns:
        (M_dict, V_dict): The computed M_k and V_k values
    """
    M, V = compute_Mk_and_Vk(degrees)
    
    # Print M_k values
    print("Mean degree of common friends to k (M_k):")
    for k, m_k in M.items():
        if not np.isnan(m_k):
            print(f"  k={k}: M_k = {m_k:.2f}")
    
    # Create figure
    plt.figure(figsize=(7, 5))
    
    # Plot histogram
    degrees_to_plot = degrees[degrees >= degree_cutoff]
    plt.hist(
        degrees_to_plot,
        bins=range(degree_cutoff, max(degrees), bin_width),
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
        # Draw small tick marks for high-degree nodes (degree > 200)
        high_degree_nodes = degrees[(degrees > 200) & (degrees < x_max)]
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
            
            # Draw line and label (M_0 gets extra height to stand out)
            if k == 0:
                plt.vlines(m_k, 0, height_M0_label, color='red', alpha=0.5)
                plt.text(m_k, height_M0_label, label, color='red',
                         ha='center', va='bottom', fontsize=14)
            else:
                plt.vlines(m_k, 0, height_Mk_labels, color='red', alpha=0.5)
                plt.text(m_k, height_Mk_labels, label, color='red',
                         ha='center', va='bottom', fontsize=14)
    
    # Format y-axis ticks for log scale
    if log_scale:
        max_power = int(np.log10(y_max)) + 1
        y_ticks = [10 ** i for i in range(max_power)]
        plt.yticks(y_ticks, [f'{tick:,.0f}' for tick in y_ticks])
    
    plt.xlabel('Degree')
    plt.ylabel('Count')
    
    return M, V


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    # Facebook network (New Orleans regional network)
    print("\n" + "="*60)
    print("NEW ORLEANS FACEBOOK NETWORK")
    print("="*60)
    
    plot_degree_distribution_with_Mk(
        NOLA.ds,
        log_scale=True,
        bin_width=1
    )
    plt.savefig(f'figures/{NOLA.name}.degree-withM.png', dpi=150, bbox_inches='tight')
    print(f"Saved: figures/{NOLA.name}.degree-withM.png")
    
    # Citation network
    print("\n" + "="*60)
    print("CITATION NETWORK")
    print("="*60)
    
    plot_degree_distribution_with_Mk(
        citations.ds,
        log_scale=True,
        show_Mk_up_to=range(0, 6),
        bin_width=1,
        m0_label_lift=0.1
    )
    plt.xlabel('In-degree (number of citations)')
    plt.savefig(f'figures/{citations.name}.degree-withM.png', dpi=150, bbox_inches='tight')
    print(f"Saved: figures/{citations.name}.degree-withM.png")
