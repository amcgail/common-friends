from common import *

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
    outf.write(
      "Mean degree of common friends to k (M_k); M_0 = population mean degree.\n"
    )
    prec = 2
    w = max(8, prec + 5)
    ks = [k for k in sorted(M.keys()) if not np.isnan(M[k])]
    m0 = M.get(0, np.nan)
    if not ks:
      outf.write("(no finite M_k)\n\n")
    else:
      outf.write(
        f"{'k':>3}  {'M_k':>{w}}  {'M_k/M_{k-1}':>{w}}  {'M_k/M_0':>{w}}\n"
        f"{'---':>3}  {'---':>{w}}  {'---':>{w}}  {'---':>{w}}\n"
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
        if np.isnan(m0) or m0 == 0:
          r0 = "—"
        else:
          r0 = f"{mk / m0:.{prec}f}"
        outf.write(f"{k:3d}  {mk_s:>{w}}  {r_prev:>{w}}  {r0:>{w}}\n")
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