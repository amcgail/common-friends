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