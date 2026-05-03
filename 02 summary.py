from common import *

_DASH = "—"


def _fmt(x, prec=2):
  """Number → string with NaN/None → '—'."""
  if x is None or (isinstance(x, float) and np.isnan(x)):
    return _DASH
  return f"{x:.{prec}f}"


def _ratio(a, b, prec=2):
  """a/b → string; '—' if either is NaN/None or b == 0."""
  if a is None or b is None:
    return _DASH
  if isinstance(a, float) and np.isnan(a):
    return _DASH
  if isinstance(b, float) and np.isnan(b):
    return _DASH
  if b == 0:
    return _DASH
  return f"{a / b:.{prec}f}"


def summ(net, tail_rows_out=None):
  if net is None:
    return
  name = net.name
  nbds = net.nbds
  nodes = net.nodes
  es = net.es
  ds = net.ds
  directed = net.directed

  # Cumulative plot: for directed citations, order cited works by in-degree (ds).
  # Previously this used out-degree of citing papers only (~65k nodes).
  if directed:
    sdeg = np.sort(ds)
  else:
    sdeg = np.array(sorted([len(nbds[n]) for n in nodes]))

  csum = np.cumsum(sdeg)
  csum = csum / csum[-1]

  plt.figure()
  plt.plot(csum)
  pct = np.arange(0, 1.1, 0.1)
  plt.yticks(pct, ["%d%%" % (p * 100) for p in pct])
  plt.ylabel("Percent of ties accounted for")
  if directed:
    plt.xlabel("Cited works (ordered by in-degree)")
  else:
    plt.xlabel("Nodes (ordered by degree)")
  plt.savefig(f"figures/{name}.degree.png")


  avg_deg = ds.mean()
  med_deg = np.median(ds)
  max_deg = ds.max()
  min_deg = ds.min()
  nmin_deg = (ds == min_deg).sum()
  var_deg = np.var(ds)

  with open(f"tables/{name}.summary.txt", 'w') as outf:
    if directed:
      out_ds = np.array([len(nbds[n]) for n in nodes])
      out_avg = out_ds.mean()
      out_var = np.var(out_ds)
      out_med = np.median(out_ds)
      n_edges = len(es)
      outf.write(
        "=== Citing papers (sources; out-degree = distinct references listed) ===\n"
        f"Count: {len(nbds)}\n"
        f"Directed citation edges: {n_edges}\n"
        f"Out-degree mean: {out_avg:0.1f}\n"
        f"Out-degree median: {out_med:0.1f}\n"
        f"Out-degree variance: {out_var:0.1f}\n"
        f"Out-degree max: {out_ds.max()}\n\n"
        "=== Cited works (targets; in-degree = times listed as a reference) ===\n"
        f"Count: {len(ds)}\n"
        f"In-degree mean: {avg_deg:0.1f}\n"
        f"In-degree median: {med_deg:0.1f}\n"
        f"In-degree variance: {var_deg:0.1f}\n"
        f"In-degree max: {max_deg}\n"
        f"In-degree min: {min_deg}\n"
        f"Cited works at min in-degree: {nmin_deg}\n"
        f"Weighted paradox analogue (in-degree): {avg_deg + var_deg / avg_deg:0.1f}\n\n"
      )
    else:
      nn = len(nbds)
      n_edges = sum(ds) / 2
      outf.write(
        "=== Nodes (undirected) ===\n"
        f"Count: {nn}\n"
        f"Edges (undirected): {n_edges:0.0f}\n"
        f"Average degree: {avg_deg:0.1f}\n"
        f"Average degree of friends: {avg_deg + var_deg / avg_deg:0.1f}\n"
        f"Median degree: {med_deg:0.1f}\n"
        f"Variance in degree: {var_deg:0.1f}\n"
        f"Max degree: {max_deg}\n"
        f"Min degree: {min_deg}\n"
        f"Nodes with min degree: {nmin_deg}\n\n"
      )

    M, _ = compute_Mk_and_Vk(ds)
    # N for Appendix E: actors whose degrees drive M_k (cited works if directed, else nodes).
    N_pop = len(ds) if directed else len(nbds)
    s2_bar = s_bar_k(2, N_pop, M)
    outf.write(
      "Mean degree of common friends to k (M_k); M_0 = population mean degree; "
      f"d_max = {max_deg:g} (maximum {('in-degree' if directed else 'degree')}).\n"
      "F(M_k): empirical CDF of degrees at M_k — F(M_k) = 100 × |{d in ds : d ≤ M_k}| / |ds|, "
      "same degree multiset ds as for M_k above.\n"
    )
    if not np.isnan(s2_bar):
      outf.write(
        f"Appendix E — approx. mean sample size to see a common friend to 2: "
        f"s̄_2 ≈ {s2_bar:.2f}  (N = {N_pop:,}, M = M_0 and V from same degree multiset).\n"
      )
    else:
      outf.write(
        "Appendix E — s̄_2: not defined (need N ≥ 2 and M^2 + V − M > 0).\n"
      )
    m2 = M.get(2, np.nan)
    if not np.isnan(m2) and len(ds) > 0:
      f_m2 = 100.0 * float((ds <= m2).mean())
      outf.write(f"Headline — F(M_2) = {f_m2:.2f}% of the population has degree ≤ M_2.\n")
    outf.write("\n")

    prec = 2
    w, w_f = 8, 10
    ks = [k for k in sorted(M.keys()) if not np.isnan(M[k])]
    if not ks:
      outf.write("(no finite M_k)\n\n")
    else:
      header = f"{'k':>3}  {'M_k':>{w}}  {'M_k/M_{k-1}':>{w}}  {'F(M_k)%':>{w_f}}  {'M_k/d_max':>{w}}\n"
      outf.write(header)
      outf.write(f"{'---':>3}  {'---':>{w}}  {'---':>{w}}  {'---':>{w_f}}  {'---':>{w}}\n")
      for k in ks:
        mk = M[k]
        f_pct = 100.0 * float((ds <= mk).mean()) if len(ds) > 0 else None
        outf.write(
          f"{k:3d}  {_fmt(mk, prec):>{w}}  {_ratio(mk, M.get(k - 1) if k > 0 else None, prec):>{w}}  "
          f"{_fmt(f_pct, prec):>{w_f}}  {_ratio(mk, max_deg, prec):>{w}}\n"
        )
      outf.write("\n")

    if len(ds) > 0:
      tail = P_k_tail_by_quantiles(ds, quantiles=(0.90, 0.95, 0.99), k_max=6)
      if tail_rows_out is not None:
        tail_rows_out.append((name, tail))
      pcts = sorted(tail.keys())
      outf.write(
        "Weighted tail exceedance P_k(d > d*), same F_{d,k} weights as M_k:\n"
        "  P_k(d > d*) = sum_{d_i > d*} F_{i,k} / sum_i F_{i,k}.\n"
        "Thresholds d* = quantile(ds, q) for q in {0.90, 0.95, 0.99} on the same multiset ds.\n\n"
      )
      for pct in pcts:
        outf.write(f"d_p{pct} = {tail[pct]['d_star']:.6g}\n")
      outf.write("\n")
      w_p = 10
      outf.write(
        f"{'k':>3}  " + "  ".join(f"{f'P_k>d{pct}':>{w_p}}" for pct in pcts) + "\n"
        f"{'---':>3}  " + "  ".join(f"{'---':>{w_p}}" for _ in pcts) + "\n"
      )
      for k in range(7):
        row_cells = [_fmt(tail[pct]["p"][k] if k < len(tail[pct]["p"]) else None, 4)
                     for pct in pcts]
        outf.write(f"{k:3d}  " + "  ".join(f"{c:>{w_p}}" for c in row_cells) + "\n")
      outf.write("\n")

    label = "In-degree" if directed else "degree"
    for i in range(2, 10):
      sub = ds[ds >= i]
      if len(sub) == 0:
        continue
      outf.write(f"Average {label} of those with {label} >= {i}: {sub.mean():0.1f}\n")

def write_P_k_tail_cross_network(tail_rows, ks=(0, 1, 2, 3), pcts=(90, 95, 99)):
  """Appendix-style table: networks × k, one block per quantile."""
  if not tail_rows:
    return
  name_w = max(len(nm) for nm, _ in tail_rows)
  col_w = 10
  lines = [
    "P_k(d > d*) with falling-factorial weights F_{d,k} (same as M_k); "
    "d* = quantile(ds, q) on the degree multiset ds used throughout the summaries.\n"
  ]
  for pct in pcts:
    lines.append(f"=== q = 0.{pct} (d_p{pct}) — P_k for k = {ks[0]}..{ks[-1]} ===\n")
    lines.append(f"{'network':<{name_w}}  " + "  ".join(f"{f'k={k}':>{col_w}}" for k in ks) + "\n")
    lines.append(f"{'-' * name_w}  " + "  ".join(f"{'-' * col_w}" for _ in ks) + "\n")
    for name, tail in tail_rows:
      ps = tail.get(pct, {}).get("p", [])
      cells = [_fmt(ps[k] if k < len(ps) else None, 4) for k in ks]
      lines.append(f"{name:<{name_w}}  " + "  ".join(f"{c:>{col_w}}" for c in cells) + "\n")
    lines.append("\n")
  Path("tables").mkdir(exist_ok=True)
  with open("tables/P_k_tail_exceedance.txt", "w") as agg:
    agg.writelines(lines)


# MAIN EXECUTION!
_tail_aggregate = []
for net in NETWORKS:
  if net is None:
    continue
  summ(net, tail_rows_out=_tail_aggregate)
write_P_k_tail_cross_network(_tail_aggregate)