"""
Build networks/highschool.txt from SocioPatterns proximity_net export.

Each raw line is one proximity record for a pupil pair (the same pair can
appear on many lines). Per-pair interaction count = number of such lines.
Let m be the median of those counts over pairs with at least one line.
Emit an undirected edge iff count > m.

Run from repo root (default paths are next to this file):
  python networks/highschool.py

Optional: PROXIMITY_CSV=path/to.csv OUT=path/to.txt python networks/...
"""
from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SRC = Path(os.environ.get("PROXIMITY_CSV", ROOT / "HighSchool2013_proximity_net.csv"))
OUT = Path(os.environ.get("OUT", ROOT / "highschool.txt"))


def main() -> None:
    if not SRC.exists():
        raise SystemExit(f"Missing input: {SRC}")

    pair_counts: dict[tuple[int, int], int] = defaultdict(int)
    with SRC.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                a, b = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            if a == b:
                continue
            if a > b:
                a, b = b, a
            pair_counts[(a, b)] += 1

    counts = np.array(list(pair_counts.values()), dtype=np.int64)
    med = float(np.median(counts))
    keep = [(i, j) for (i, j), c in pair_counts.items() if c > med]
    keep.sort()

    with OUT.open("w") as w:
        for i, j in keep:
            w.write(f"{i} {j}\n")

    print(
        f"wrote {OUT}: {len(keep)} edges "
        f"(count > median={med:g}; {len(pair_counts)} pairs with ≥1 row, max {counts.max()})"
    )


if __name__ == "__main__":
    main()
