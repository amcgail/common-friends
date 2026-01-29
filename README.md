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
| `01 create citation network.py` | *Optional.* Builds `networks/citations.txt` from Web of Science export files (see below). |
| `02 summary.py` | Network summary statistics and degree-distribution figures; computes \(M_0, M_1, \ldots, M_5\) from the degree distribution. |
| `03 count common friends.py` | Simulations: sample until common friend to \(k\) for \(k = 2,\ldots,6\); 300 runs; writes tables with sample size and degree (or indegree) statistics. |
| `networks/` | Edge-list inputs: `neworleans.txt`, `citations.txt`. |
| `tables/` | Summary and simulation output (e.g. `*.summary.txt`, `*.sim.txt`). |
| `figures/` | Degree-distribution and related plots. |

## Requirements

- Python 3.7+
- Dependencies: `numpy`, `matplotlib`, `tqdm`; optional `python-dotenv` for `01 create citation network.py`.

Install with:

```bash
pip install -r requirements.txt
```

## How to run (replication)

1. **Networks**: Ensure `networks/neworleans.txt` and `networks/citations.txt` exist (see *Data* below).  
   If you already have these files, skip step 2.

2. **Citation network (optional)**  
   Only needed if you are building the citation network from Web of Science exports.  
   Set environment variable `BASE_DIR` to the folder containing the exported `.txt` files, then run:
   ```bash
   python "01 create citation network.py"
   ```
   This writes `networks/citations.txt`. If you are only replicating from existing edge lists, use the provided `citations.txt` and do not run this script.

3. **Summary statistics and \(M_k\)**  
   ```bash
   python "02 summary.py"
   ```  
   Writes `tables/<network>.summary.txt` and `figures/<network>.degree.png` for each network. The summary includes \(M_0\) (mean degree), \(M_1 = M_0 + V_0/M_0\) (mean degree of friends), and \(M_2, \ldots, M_5\) from the recurrence in Equation 1.

4. **Common-friend simulations**  
   ```bash
   python "03 count common friends.py"
   ```  
   For each network, samples nodes until observing a common friend to \(k\) for \(k = 2, 3, 4, 5, 6\), over 300 runs. Writes `tables/<network>.sim.txt` with mean and 5%–95% quantiles for sample size and degree (Facebook) or indegree (citations), matching the design described in the paper (Methods: Simulations).

## Data

- **New Orleans Facebook**: From Viswanath et al. (2009). The edge list `networks/neworleans.txt` should contain one edge per line: `node_id_i node_id_j` (integer IDs). The code symmetrizes this network so that “friends” are undirected.
- **Citations**: From Web of Science (sociology journals, 2010–2019). The edge list `networks/citations.txt` should be one citation per line: `citing_paper_id cited_work_id`. The code treats this as directed (outdegree = citations made, indegree = citations received); the paper’s “common friends” in this setting are *sources cited in common*, i.e. high-*indegree* nodes.
