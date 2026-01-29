from common import *

def sample_until_common(net, N=2):
  common = defaultdict(int)
  chosen = set()
  buffer = None
  while not len(common) or max(common.values()) < N:
    if buffer is None or not len(buffer):
      candidates = [n for n in net.nodes if n not in chosen]
      buffer = list(np.random.choice(candidates, 200, replace=False))

    new_node = buffer.pop()

    chosen.add(new_node)

    for nb in net.nbds[new_node]:
      common[nb] += 1

  max_degree = max(common.values())
  candidates = [n for n in common if common[n] == max_degree]
  result = np.random.choice(candidates)
  return len(chosen), result

def sim_study(net):
  print('Starting simulation for', net.name)
  Ns = [2,3,4,5,6]
  
  d = [len(n) for n in net.nbds.values()]
  d = np.array(d)
  max_degree = d.max()
  top_degrees = d.sort()[-5:]

  stats_wait = defaultdict(list)
  for N in Ns:
    print('working on', N)
    samples = []
    samples_ = []
    for i in range(300):
      d_, ni = sample_until_common(net, N=N)
      d = len(net.nbds[ni])

      samples.append(d)
      samples_.append(d_)

    stats_wait["avg_deg"].append(np.mean(samples))
    stats_wait["max_deg"].append(np.max(samples))
    stats_wait["min_deg"].append(np.min(samples))
    stats_wait["med_deg"].append(np.median(samples))
    stats_wait["avg_samp"].append(np.mean(samples_))
    stats_wait["5pct_samp"].append(np.percentile(samples_, 5))
    stats_wait["95pct_samp"].append(np.percentile(samples_, 95))
    stats_wait["5pct_deg"].append(np.percentile(samples, 5))
    stats_wait["95pct_deg"].append(np.percentile(samples, 95))
    stats_wait["got_max"].append(np.mean([s==max_degree for s in samples]))
    stats_wait["got_top5"].append(np.mean([s in top_degrees for s in samples]))

  with open(f'tables/{net.name}.sim.txt', 'w') as outf:
    for ni, N in enumerate(Ns):
      outf.write(f'------ {N} ------\n')
      for stat in stats_wait:
          outf.write(f'{stat}: {stats_wait[stat][ni]}\n')
      outf.write(f'avg_deg percentile: {1-(ds >= stats_wait["avg_deg"][ni]).sum() / len(ds):0.3%}\n')

sim_study(NOLA)
sim_study(citations)