from common import *


def approx_s2_bar(N, M, V):
  """Appendix E: mean sample size to observe a common friend to 2, ~ 1/2 + sqrt(1/4 + 2(N-1)/(M^2+V-M))."""
  if N is None or N < 2 or np.isnan(M) or M <= 0:
    return np.nan
  denom = M * M + V - M
  if denom <= 0 or np.isnan(denom):
    return np.nan
  return 0.5 + np.sqrt(0.25 + 2 * (N - 1) / denom)


def summ(net):
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
    s2_bar = approx_s2_bar(N_pop, float(avg_deg), float(var_deg))
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
    w = max(8, prec + 5)
    w_f = max(10, prec + 6)
    ks = [k for k in sorted(M.keys()) if not np.isnan(M[k])]
    if not ks:
      outf.write("(no finite M_k)\n\n")
    else:
      outf.write(
        f"{'k':>3}  {'M_k':>{w}}  {'M_k/M_{k-1}':>{w}}  {'F(M_k)%':>{w_f}}  {'M_k/d_max':>{w}}\n"
        f"{'---':>3}  {'---':>{w}}  {'---':>{w}}  {'---':>{w_f}}  {'---':>{w}}\n"
      )
      for k in ks:
        mk = M[k]
        mk_s = f"{mk:.{prec}f}"
        if k == 0:
          r_prev = "—"
        else:
          m_prev = M.get(k - 1, np.nan)
          if np.isnan(m_prev) or m_prev == 0:
            r_prev = "—"
          else:
            r_prev = f"{mk / m_prev:.{prec}f}"
        if len(ds) == 0:
          f_pct = "—"
        else:
          f_pct = f"{100.0 * float((ds <= mk).mean()):.{prec}f}"
        if max_deg == 0:
          rdmax = "—"
        else:
          rdmax = f"{mk / max_deg:.{prec}f}"
        outf.write(f"{k:3d}  {mk_s:>{w}}  {r_prev:>{w}}  {f_pct:>{w_f}}  {rdmax:>{w}}\n")
      outf.write("\n")

    label = "In-degree" if directed else "degree"
    for i in range(2, 10):
      sub = ds[ds >= i]
      if len(sub) == 0:
        continue
      outf.write(f"Average {label} of those with {label} >= {i}: {sub.mean():0.1f}\n")

# MAIN EXECUTION!
for net in NETWORKS:
    summ(net)