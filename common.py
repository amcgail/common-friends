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

# --- M_k / V_k and tail probabilities from a degree multiset (Appendix B). ---


def _F_weights(degrees, k_max):
  """
  Yield (k, w) for k = 0..k_max where w[i] = F_{d_i,k} = d_i (d_i-1) ... (d_i-k+1).
  k = 0 → w = ones; subsequent k builds on the previous w (one O(N) multiply per k).
  Same falling-factorial weights used to define M_k.
  """
  d = np.asarray(degrees, dtype=np.float64)
  w = np.ones_like(d)
  yield 0, w
  for k in range(1, k_max + 1):
    w = w * (d - (k - 1))
    np.maximum(w, 0.0, out=w)
    yield k, w


def compute_Mk_and_Vk(degrees, k_max=24):
  """Weighted mean/variance of degree under F_{d,k} (Appendix B, Equations 7–8)."""
  if len(degrees) == 0:
    return {0: np.nan}, {0: np.nan}
  d = np.asarray(degrees, dtype=np.float64)
  M, V = {}, {}
  for k, w in _F_weights(d, k_max):
    s = float(w.sum())
    if s <= 0:
      M[k] = V[k] = np.nan
      continue
    mk = float((w * d).sum() / s)
    M[k] = mk
    V[k] = float((w * (d - mk) ** 2).sum() / s)
  return M, V


def E_k_expected(s, k, N, M_dict):
  """
  Appendix E: expected count of size-k subsets of an s-sample with a common
  friend in the population, via the recurrence
      E_1 = s * M_0,    E_{k+1} = E_k * (s - k) / (N - k) * (M_k - k) / (k + 1).
  Uses M_dict (e.g. from compute_Mk_and_Vk). Returns 0 if any factor is non-positive
  (e.g. s < k); NaN if the needed M_k is missing.
  """
  if k < 1 or s <= 0 or N is None or N < k:
    return 0.0
  M0 = M_dict.get(0, np.nan)
  if np.isnan(M0):
    return np.nan
  val = s * M0
  for j in range(1, k):
    Mj = M_dict.get(j, np.nan)
    if np.isnan(Mj):
      return np.nan
    factor = (s - j) * (Mj - j) / ((N - j) * (j + 1))
    if factor < 0:
      return 0.0
    val *= factor
  return val


def s_bar_k(k, N, M_dict, target=1.0):
  """Smallest real s with E_k(s) >= target (Appendix E inversion). NaN if unattainable."""
  if k < 1 or N is None or N < 1 or not M_dict:
    return np.nan
  if E_k_expected(float(N), k, N, M_dict) < target:
    return np.nan
  lo, hi = 0.0, float(N)
  for _ in range(200):
    mid = 0.5 * (lo + hi)
    if E_k_expected(mid, k, N, M_dict) < target:
      lo = mid
    else:
      hi = mid
    if hi - lo < 1e-12 * max(1.0, hi):
      break
  return 0.5 * (lo + hi)


def P_k_tail_by_quantiles(degrees, quantiles=(0.90, 0.95, 0.99), k_max=6):
  """
  P_k(d > d*) = (sum_i F_{i,k} 1[d_i > d*]) / sum_i F_{i,k},
  with d* = quantile(degrees, q) on the same multiset. Returns
  {pct: {"d_star": d_star, "p": [P_0, ..., P_{k_max}]}}.
  """
  d = np.asarray(degrees, dtype=np.float64)
  if len(d) == 0:
    return {}
  thresholds = {int(round(100 * q)): float(np.quantile(d, q)) for q in quantiles}
  out = {pct: {"d_star": dstar, "p": []} for pct, dstar in thresholds.items()}
  for _, w in _F_weights(d, k_max):
    s = float(w.sum())
    for pct, dstar in thresholds.items():
      out[pct]["p"].append(np.nan if s <= 0 else float(w[d > dstar].sum() / s))
  return out


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