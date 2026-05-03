# Common Friends: Replication Code

Replication code for **"The Power of Common Friends"** (McGail & Feld).  
This repository reproduces the empirical illustrations and summary statistics reported in the paper.

## What the paper shows

- **Friendship paradox**: On average, your friends have more friends than you do; mean degree of friends is \(M_f = M + V/M\) (population mean \(M\), variance \(V\)).
- **Power of common friends**: A *common friend to k* is someone who appears in the friend lists of \(k\) different people. The mean degree of common friends to \(k\) increases with \(k\); the paper gives the recurrence \(M_{k+1} = M_k + V_k/(M_k - k)\) (Equation 1 in the main text).
- **Empirical illustrations**:  
  - **Facebook (New Orleans)**: Sampling people until we see a common friend to 2, 3, … shows that those common friends have very high degree (e.g. common friend to 3 has mean degree &gt; 99.9th percentile).  
  - **Citations**: Sampling papers until we see a source cited by 2, 3, … papers shows that those sources have very high citation counts (indegree).

## Repository layout

| File / folder | Purpose |
|---------------|--------|
| `common.py` | Shared `Network` class, degree/neighborhood handling, and computation of \(M_k\) from the degree distribution (Equation 1). |
| `02 summary.py` | Network summary statistics and degree-distribution figures; computes \(M_0, M_1, \ldots, M_5\) from the degree distribution. |
| `03 count common friends.py` | Simulations: sample until common friend to \(k\) for \(k = 2,\ldots,6\); 300 runs; writes tables with sample size and degree (or indegree) statistics. |
| `networks/` | `*.txt` edge lists; same-folder `*.py` files generate the matching `.txt` where applicable. |
| `tables/` | Summary and simulation output (e.g. `*.summary.txt`, `*.sim.txt`). |
| `figures/` | Degree-distribution and related plots. |

## Requirements

- Python 3.7+
- Dependencies: `numpy`, `matplotlib`, `tqdm`; optional `python-dotenv` for `networks/citations.py` (WoS).

Install with:

```bash
pip install -r requirements.txt
```

## How to run (replication)

1. **Networks**: Ensure the `*.txt` edge lists you need exist under `networks/` (see *Data* below); missing files are skipped in `common.py`.

2. **Summary statistics and \(M_k\)**  
   ```bash
   python "02 summary.py"
   ```  
   Writes `tables/<network>.summary.txt` and `figures/<network>.degree.png` for each network. The summary includes \(M_0\) (mean degree), \(M_1 = M_0 + V_0/M_0\) (mean degree of friends), and \(M_2, \ldots, M_5\) from the recurrence in Equation 1.

3. **Common-friend simulations**  
   ```bash
   python "03 count common friends.py"
   ```  
   For each network, samples nodes until observing a common friend to \(k\) for \(k = 2, 3, 4, 5, 6\), over 300 runs. Writes `tables/<network>.sim.txt` with mean and 5%–95% quantiles for sample size and degree (Facebook) or indegree (citations), matching the design described in the paper (Methods: Simulations).

## Data

Edge lists are **one line per edge**, **space-separated** integer node IDs. Extra columns (if any) are ignored by the loader.

| File | Directed? | In `common.py` | Source |
|------|-----------|----------------|--------|
| `neworleans.txt` | No (symmetrized) | Yes | Facebook New Orleans network, Viswanath et al. (2009). |
| `citations.txt` | Yes | Yes | Web of Science, sociology journals 2010–2019; `citing_id cited_id`. |
| `indian_microfinance.txt` | No (symmetrized) | Yes | Banerjee et al. (2013), [doi:10.7910/DVN/U3BIHX](https://doi.org/10.7910/DVN/U3BIHX); village 60, all-village undirected layer. |
| `power.txt` | No (symmetrized) | Yes | Western US power grid, Watts & Strogatz (1998); from SuiteSparse “Newman/power” Matrix Market (`power.mtx`). Nodes **0 … 4940** (MTX 1-based converted). |
| `highschool.txt` | No (symmetrized) | Yes | SocioPatterns Marseille high school 2013 proximity (`HighSchool2013_proximity_net.csv`). Each raw row is one proximity record for a pair (the same pair can appear many times). **Per-pair interaction count** = number of such rows. An undirected edge is kept iff that count is **strictly greater than the median** of those counts over all pairs with at least one row. [Dataset page](https://sociopatterns.org/datasets/high-school-contact-and-friendship-networks/). |
| `celegans_chemical.txt` | Yes | Yes | *C. elegans* hermaphrodite **chemical** connectome, Cook et al. (2019), corrected matrices via [Netzschleuder `celegans_2019`](https://networks.skewed.de/net/celegans_2019) / [WormWiring](https://wormwiring.org/pages/adjacency.html); edges from `hermaphrodite_chemical_corrected.csv/edges.csv` (pre→post). |

- **New Orleans Facebook**: `neworleans.txt` may use tabs; `node_id_i node_id_j` (1-based in the original release). Symmetrized so neighborhoods are undirected.
- **Citations**: `citing_paper_id cited_work_id`. Directed: indegree = times cited; simulations use indegree for the “common friend” target.
- **Indian microfinance**: One undirected edge per line `i j` (0-based indices). Symmetrized like New Orleans.
