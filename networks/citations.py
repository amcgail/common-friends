"""
Build networks/citations.txt from Web of Science plain-text exports.

Run from anywhere (repo root recommended):
  BASE_DIR=/path/to/wos/exports python networks/citations.py

If you are replicating from the bundled citations.txt, you do not need this.
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()

BASE_DIR = os.environ.get("BASE_DIR")
if not BASE_DIR:
    raise SystemExit("Set BASE_DIR to the folder containing WoS exported .txt files.")
BASE_DIR = Path(BASE_DIR)
print("base dir:", BASE_DIR)

fns = list(BASE_DIR.glob("**/*.txt"))
print(len(fns), "files found")

work_range = 2010, 2019  # citing population


def normalize_c(x):
    a, *other = x.split(", ")
    aparts = a.split()
    if len(aparts) > 1:
        l, f, *_ = aparts
        f = f[0]
        a = f"{l} {f}"

    if not (re.search(r"\b(v|p)[0-9]+\b", x) or "doi" in x):
        try:
            int(other[0])
            other = other[1:]
        except Exception:
            pass

    return ", ".join([a] + other)


cdict = {}
skipped_no_date = 0
for fn in tqdm(fns):
    with open(fn) as f:
        for li, l in enumerate(f):
            if li == 0:
                head = l.split()
                continue
            parts = l.split("\t")
            py, doi = parts[head.index("PY")], parts[head.index("DI")]
            if py == "":
                if doi != "":
                    py = doi.split(".")
                    py = [x for x in py if x.isdigit() and 1900 < int(x) < 2020]
                    if len(py):
                        py = int(py[0])
                    else:
                        py = ""
                else:
                    py = ""
            else:
                py = int(py)

            if py == "":
                if parts[66] != "":
                    py = int(parts[66].split("-")[0])

            if py == "":
                skipped_no_date += 1
                continue

            if py < work_range[0] - 1 or py > work_range[1] + 1:
                break

            if py < work_range[0] or py > work_range[1]:
                continue

            citations = parts[29].split(";")
            citations = [x.strip().lower() for x in citations if x.strip() != ""]
            citations = [normalize_c(x) for x in citations]
            citations = [x for x in citations if x != "[anonymous], communication"]
            paper_id = parts[61]

            cdict[paper_id] = citations

for x in cdict:
    cdict[x] = list(set(cdict[x]))

index = {}


def i(name):
    if name not in index:
        index[name] = len(index)
    return index[name]


out_path = REPO / "networks" / "citations.txt"
with out_path.open("w") as outf:
    for f, cits in cdict.items():
        for t in cits:
            outf.write(f"{i(f)} {i(t)}\n")

print("wrote", out_path, "skipped_no_date:", skipped_no_date)
