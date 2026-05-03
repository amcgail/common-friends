from pathlib import Path
from collections import defaultdict, Counter
import re
import os
from tqdm import tqdm
from dotenv import load_dotenv
from matplotlib import pyplot as plt
import numpy as np

Path('figures').mkdir(exist_ok=True)
Path('tables').mkdir(exist_ok=True)

class Network:
  def __init__(self, fn, symmetrize=False):
    self.name = Path(fn).stem
    self.loaded = False
    if not Path(fn).exists():
      print('Warning:', fn, 'doesn\'t exist')
      return

    # loop through edgelist, gathering neighborhoods as we go
    self.nbds = defaultdict(set)
    self.es = []
    with open(fn) as f:
        for l in f:
            ai,bi,*_ = l.split()
            ai,bi = int(ai),int(bi)
            
            if symmetrize: self.nbds[bi].add(ai) # for the new orleans network
            self.nbds[ai].add(bi)
            self.es.append((ai,bi))

    self.directed = not symmetrize
    self.nodes = list(self.nbds.keys())
    # Multiset of degrees: undirected = incident edges per node; directed = in-degree per cited work
    self.ddict = Counter([y for n in self.nbds.values() for y in n])
    self.ds = np.array(list(self.ddict.values()))
    self.loaded = True


def _network_or_none(path, symmetrize=False):
    net = Network(path, symmetrize=symmetrize)
    return net if net.loaded else None


# DEFINE THE NETWORKS (omit entries when the edge list is missing)
NOLA = _network_or_none('networks/neworleans.txt', symmetrize=True)
citations = _network_or_none('networks/citations.txt', symmetrize=False)
indian_microfinance = _network_or_none('networks/indian_microfinance.txt', symmetrize=True)
power = _network_or_none('networks/power.txt', symmetrize=True)
highschool = _network_or_none('networks/highschool.txt', symmetrize=True)
celegans_chemical = _network_or_none('networks/celegans_chemical.txt', symmetrize=False)

NETWORKS = [
  NOLA,
  citations,
  indian_microfinance,
  power,
  highschool,
  celegans_chemical
]

# --- M_k / V_k from degree multiset (Appendix B, Equations 7–8) ---

_falling_factorial_cache = {}


def falling_factorial(d, k):
  """F_{d,k} = d(d-1)...(d-k+1); zero if d < k (handled by product)."""
  if (d, k) in _falling_factorial_cache:
    return _falling_factorial_cache[d, k]
  result = 1
  for i in range(k):
    result *= (d - i)
  _falling_factorial_cache[d, k] = result
  return result


def falling_factorial_array(max_degree, k):
  """F_{d,k} for d = k..max_degree; index d - k into the returned array."""
  F_k_k = falling_factorial(k, k)
  weights = [F_k_k]
  for d in range(k + 1, max_degree + 1):
    weights.append(weights[-1] * d / (d - k))
  return np.array(weights)


def compute_Mk_and_Vk(degrees, max_k=25):
  """Weighted mean/variance of degree for common-friends-to-k weights F_{d,k}."""
  if len(degrees) == 0:
    return {0: np.nan}, {0: np.nan}
  M = {0: degrees.mean()}
  V = {0: degrees.var()}
  max_d = int(max(degrees))
  for k in range(1, max_k):
    weights_by_degree = falling_factorial_array(max_d, k)
    eligible_degrees = degrees[degrees >= k]
    if len(eligible_degrees) == 0:
      M[k] = np.nan
      V[k] = np.nan
      continue
    node_weights = weights_by_degree[eligible_degrees - k]
    M[k] = (node_weights * eligible_degrees).sum() / node_weights.sum()
    V[k] = (node_weights * (eligible_degrees - M[k]) ** 2).sum() / node_weights.sum()
  return M, V


def format_mk_lines(M, indent="", precision=2):
  lines = []
  for k in sorted(M.keys()):
    m_k = M[k]
    if np.isnan(m_k):
      continue
    lines.append(f"{indent}k={k}: M_k = {m_k:.{precision}f}")
  return lines


def print_mk(M, title="Mean degree of common friends to k (M_k):"):
  print(title)
  for line in format_mk_lines(M, indent="  ", precision=2):
    print(line)