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
            
            if symmetrize: self.nbds[ai].add(bi) # for the new orleans network
            self.nbds[bi].add(ai)
            self.es.append((ai,bi))

    self.nodes = list(self.nbds.keys())
    self.name = Path(fn).stem

# DEFINE THE NETWORKS
NOLA = Network('networks/neworleans.txt', symmetrize=True)
citations = Network('networks/citations.txt', symmetrize=False)