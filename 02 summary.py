from common import *

def summ(net):
  name = net.name
  nbds = net.nbds
  nodes = net.nodes
  es = net.es
  ds = net.ds

  # CUMULATIVE PLOT OF DEGREE!
  sdeg = np.array(sorted([len(nbds[n]) for n in nodes]))
  csum = np.cumsum(sdeg)
  csum = csum / csum[-1]

  plt.figure()
  plt.plot(csum)
  # percentages
  pct = np.arange(0,1.1,0.1)
  plt.yticks(pct,["%d%%" % (p*100) for p in pct]);
  plt.xlabel("Nodes")
  plt.ylabel("Percent of ties accounted for")
  plt.savefig(f"figures/{name}.degree.png")


  # BASIC STATISTICS
  nn = len(nbds)
  
  avg_deg = ds.mean()
  med_deg = np.median(ds)
  max_deg = ds.max()
  min_deg = ds.min()
  nmin_deg = (ds==min_deg).sum()
  var_deg = np.var(ds)

  max_degree = max(len(nbds[i]) for i in nbds)
  top_degrees = set(sorted([len(nbds[i]) for i in nbds], reverse=True)[:5])

  with open(f"tables/{name}.summary.txt", 'w') as outf:
    outf.write(f"""
Number of nodes: {nn}
Number of edges: {sum(ds)/2}
Average degree: {avg_deg:0.1f}
Average degree of friends: {avg_deg + var_deg/avg_deg:0.1f}
Median degree: {med_deg:0.1f}
Variance in degree: {var_deg:0.1f}
Max degree: {max_deg}
Min degree: {min_deg}
Number of nodes with min degree: {nmin_deg}
    \n""")

    for i in range(2,10):
        outf.write(f"Average degree of those with degree >={i}: {ds[ds>=i].mean():0.1f}\n")

# MAIN EXECUTION!
summ(NOLA)
summ(citations)