"""
05 tradeoff figure.py
=====================
Per network, the analytical trade-off between

  x = expected sample size  s̄_k  to observe one common friend to k   (Appendix E)
  y = P_k(d > d_pQ)  with falling-factorial weights F_{i,k}           (Appendix B)

for k = 1..15 and Q in {90, 95, 99}. Horizontal dashed line marks the random-
sampling baseline (1 − q). All quantities are deterministic functions of the
degree distribution; no Monte Carlo.

Outputs (per quantile Q in {90, 95, 99}):
  figures/tradeoff_p{Q}_frac.png — x = s̄_k / N (percent tick labels)
  figures/tradeoff_p{Q}_sbar.png — x = s̄_k (log scale)
  tables/tradeoff.tsv   (network  N  k  s_bar_k  s_bar_frac  P_k_p90  P_k_p95  P_k_p99)
"""

from common import *

from matplotlib.ticker import FuncFormatter, LogFormatterSciNotation


def _pct_label(x, _pos):
  """Format a fraction as a percent label without scientific notation (e.g. 1e-6 -> '0.0001%')."""
  if x <= 0:
    return ""
  s = f"{x * 100:.10f}".rstrip("0").rstrip(".")
  return f"{s}%"


def _trim_curve(ks, xs, ys, p_cap=0.99, k_min_keep=4):
  """Keep points up to (but not including) the first y strictly above p_cap;
  always keep at least all points with k ≤ k_min_keep (and ≥ 1 point)."""
  end = len(ys)
  for i, p in enumerate(ys):
    if p > p_cap:
      end = i
      break
  min_end = 0
  for i, k in enumerate(ks):
    if k <= k_min_keep:
      min_end = i + 1
  end = max(end, min_end, 1)
  return ks[:end], xs[:end], ys[:end]


DISPLAY = {
  "neworleans": "Facebook (New Orleans)",
  "citations": "Sociology citations",
  "indian_microfinance": "Indian microfinance",
  "power": "Western US power grid",
  "highschool": "Marseille high school",
  "celegans_chemical": "C. elegans chemical",
}

QUANTILES = (0.90, 0.95, 0.99)
K_VALUES = tuple(range(1, 16))
PCTS = tuple(int(round(100 * q)) for q in QUANTILES)

LEGEND_ORDER = (
  "citations",
  "neworleans",
  "indian_microfinance",
  "power",
  "celegans_chemical",
  "highschool",
)
# Draw C. elegans before the power grid so the power curve (markers + line) is on top where they overlap.
DRAW_ORDER = (
  "citations",
  "neworleans",
  "indian_microfinance",
  "celegans_chemical",
  "power",
  "highschool",
)
def network_data(net):
  """For one network: aligned ks, s_bar_k, s_bar_k/N, and P_k per quantile (NaN points dropped)."""
  ds = net.ds
  N, M_base = network_population_and_M_base(net)
  M, _ = compute_Mk_and_Vk(ds, k_max=max(K_VALUES))
  tail = P_k_tail_by_quantiles(ds, quantiles=QUANTILES, k_max=max(K_VALUES))
  data = {"name": net.name, "N": N, "M_base": M_base, "ks": [], "s": [], "frac": [],
          "p": {pct: [] for pct in PCTS}}
  for k in K_VALUES:
    s = s_bar_k(k, N, M, M_base=M_base)
    if np.isnan(s):
      continue
    ps = {pct: tail.get(pct, {}).get("p", [])[k] if k < len(tail.get(pct, {}).get("p", [])) else np.nan
          for pct in PCTS}
    if any(np.isnan(p) for p in ps.values()):
      continue
    data["ks"].append(k)
    data["s"].append(float(s))
    data["frac"].append(float(s) / N)
    for pct, p in ps.items():
      data["p"][pct].append(float(p))
  return data


def plot_tradeoff(all_data, pct, use_frac):
  """use_frac: True → x = s̄_k/N with percent ticks; False → x = s̄_k (log, default formatter)."""
  baseline = 1.0 - pct / 100.0
  xkey = "frac" if use_frac else "s"
  # PNAS double-column width (7 in) keeps every curve in the legend distinguishable.
  fig, ax = plt.subplots(figsize=(PNAS_WIDTHS_IN["2col"], 5.0))
  lines_by_name = {}
  for d in all_data:
    if not d["ks"]:
      continue
    ks, xs, ys = _trim_curve(d["ks"], d[xkey], d["p"][pct], p_cap=0.995)
    (line,) = ax.plot(xs, ys, linestyle="-", linewidth=1.4, color=None, zorder=1)
    lines_by_name[d["name"]] = line
    c = line.get_color()
    # Draw markers + labels in descending k so each (disk, digit) stacks: high k
    # underneath, low k on top; text sits just above its own disk (z+1).
    for j, idx in enumerate(reversed(range(len(ks)))):
      zz = 4 + 2 * j
      ax.scatter(
        [xs[idx]],
        [ys[idx]],
        s=78,
        facecolors="white",
        edgecolors=c,
        linewidths=2.0,
        zorder=zz,
        clip_on=False,
      )
      ax.annotate(
        str(ks[idx]),
        (xs[idx], ys[idx]),
        ha="center",
        va="center",
        fontsize=6.25,
        fontweight="bold",
        color="#1a1a1a",
        zorder=zz + 1,
        clip_on=False,
      )
  hline = ax.axhline(
    baseline,
    linestyle="--",
    color="gray",
    alpha=0.7,
    zorder=1.5,
    label=f"random-sample baseline ({baseline:.2f})",
  )
  ax.set_xscale("log")
  if use_frac:
    ax.xaxis.set_major_formatter(FuncFormatter(_pct_label))
    ax.set_xlabel(r"Sample fraction  $\bar s_k / N$  to find one common friend to $k$")
  else:
    ax.xaxis.set_major_formatter(LogFormatterSciNotation())
    ax.set_xlabel(r"Expected sample size  $\bar s_k$  to find one common friend to $k$")
  ax.set_ylabel(rf"$P_k(d > d_{{p{pct}}})$  (analytical, $F_{{i,k}}$ weights)")
  ax.set_ylim(0.0, 1.02)
  ax.grid(True, which="both", linestyle=":", alpha=0.4)
  leg_handles = [lines_by_name[n] for n in LEGEND_ORDER if n in lines_by_name]
  leg_labels = [DISPLAY[n] for n in LEGEND_ORDER if n in lines_by_name]
  leg_handles.append(hline)
  leg_labels.append(hline.get_label())
  leg = ax.legend(leg_handles, leg_labels, fontsize=8, loc="lower right", framealpha=0.85)
  if leg is not None:
    leg.set_zorder(10_000)
  fig.tight_layout()
  suffix = "frac" if use_frac else "sbar"
  written = save_figure(f"tradeoff_p{pct}_{suffix}", fig=fig)
  plt.close(fig)
  for path in written:
    print(f"Saved: {path}")


def main():
  Path("figures").mkdir(exist_ok=True)
  Path("tables").mkdir(exist_ok=True)

  all_data = [network_data(net) for net in NETWORKS if net is not None]
  draw_idx = {name: i for i, name in enumerate(DRAW_ORDER)}
  legend_idx = {name: i for i, name in enumerate(LEGEND_ORDER)}
  data_draw_order = sorted(all_data, key=lambda d: draw_idx.get(d["name"], len(DRAW_ORDER)))
  data_legend_order = sorted(all_data, key=lambda d: legend_idx.get(d["name"], len(LEGEND_ORDER)))

  for pct in PCTS:
    plot_tradeoff(data_draw_order, pct, use_frac=True)
    plot_tradeoff(data_draw_order, pct, use_frac=False)

  out_tsv = "tables/tradeoff.tsv"
  with open(out_tsv, "w") as f:
    f.write(
      f"# s_bar_k from Appendix E; P_k(d > d_pQ) with F_{{i,k}} weights;"
      f" k = {K_VALUES[0]}..{K_VALUES[-1]}.\n"
    )
    f.write("network\tN\tk\ts_bar_k\ts_bar_frac\t" + "\t".join(f"P_k_p{pct}" for pct in PCTS) + "\n")
    for d in data_legend_order:
      for i, k in enumerate(d["ks"]):
        ps = "\t".join(f"{d['p'][pct][i]:.6g}" for pct in PCTS)
        f.write(f"{d['name']}\t{d['N']}\t{k}\t{d['s'][i]:.6g}\t{d['frac'][i]:.6g}\t{ps}\n")
  print(f"Saved: {out_tsv}")


if __name__ == "__main__":
  main()
