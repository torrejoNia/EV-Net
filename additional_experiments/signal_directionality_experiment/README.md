# Signal directionality experiment

EV-Net propagates signals through the NicheNet ligand–signaling network
(`lr_sig_network`) and uses **all nodes** of that network as possible seed
nodes, not only ligands. One possible concern with this design is the
direction of the signal: could a random walk started from a downstream
signaling node travel *back* to ligands instead of moving downstream toward
transcriptional regulators?

These three experiments test that on the NicheNet **mouse** networks, using
the Python package [igraph](https://python.igraph.org/).

**Result in short:** ligands that are not also receptors have no incoming
edges in `lr_sig_network`. Downstream nodes can therefore never reach them,
either by any directed path or by a random walk. Using all nodes as seeds
does not reverse the direction of the signal.

## Contents

| File | Description |
|---|---|
| `experiment1_shortest_path.py` | Exp. 1: shortest paths from ligands vs. downstream nodes to transcriptional regulators |
| `experiment2_random_walk.py` | Exp. 2: random walk with restart (personalized PageRank) started from downstream nodes |
| `experiment3_reachability_to_ligands.py` | Exp. 3: can downstream nodes reach ligands by any directed path? |
| `experiment3_supportgraph.py` | Toy network showing how directed reachability works in igraph (support for Exp. 3) |
| `lr_network.parquet` | NicheNet mouse ligand → receptor network |
| `sig_network.parquet` | NicheNet mouse signaling network |
| `lr_sig_network.parquet` | Weighted ligand–signaling network used as the graph in all experiments (`from`, `to`, `weight`) |
| `figures/` | Shortest path histograms produced by Exp. 1 |

The networks were exported from the NicheNet R objects to
[Parquet](https://arrow.apache.org/docs/python/parquet.html) for use from
Python. The weighted network has 19,763 nodes and 3,865,137 directed edges.

## How to run

Tested with Python 3.10.12.

```bash
cd additional_experiments/signal_directionality_experiment
pip install -r requirements.txt
python experiment1_shortest_path.py
python experiment2_random_walk.py
python experiment3_reachability_to_ligands.py
python experiment3_supportgraph.py
```

The scripts are split into `#%%` cells, so they can also be run cell by cell
in VS Code or Spyder. Run them from this folder, because the data files are
read with relative paths.

- Experiments 1 and 3 have a `USE_UNIT_WEIGHTS` switch at the top:
  `False` uses the original NicheNet weights, `True` sets every edge weight
  to 1.
- Random samples use `random.seed(42)`.
- Each script ends with a **"How to read these metrics"** section.
- Running time on a 16-core laptop: Exp. 1 takes about 1 min per weight
  setting, Exp. 3 a few minutes, and Exp. 2 (300 random walks) up to about
  an hour.

## Node categories

The network is highly interconnected and many genes have more than one
role, so the categories differ slightly between experiments:

| Category | Exp. 1 | Exp. 2 and 3 |
|---|---|---|
| Ligands | `lr_network$from` (1,287; 1,000 sampled) | `lr_network$from` minus receptors (1,005) |
| Receptors | – | `lr_network$to` (1,084) |
| Downstream nodes | `sig_network$from` minus ligands (18,405; 1,000 sampled) | nodes of `sig_network` (`from` or `to`) absent from `lr_network` (17,674) |
| Transcriptional regulators | `sig_network$to` (18,729) | – |

282 genes act both as ligand and as receptor. In experiments 2 and 3 they
are counted as receptors only.

## Results

### 1. Shortest path analysis

Average shortest path length to transcriptional regulators:

| Weights | Ligands → TRs | Downstream → TRs |
|---|---|---|
| Original | 0.22 | 0.21 |
| All = 1 | 2.40 | 2.43 |

Ligands and downstream nodes are equally far from the transcriptional
regulators: almost every node reaches almost every regulator in 2–3 steps
(about 99.5% of pairs are reachable). This reflects the small-world
structure of the network, so shortest paths alone cannot show directionality.

Note: igraph treats weights as costs, while in NicheNet a higher weight
means a *more* reliable interaction. The weighted values are therefore
useful only for comparing the two groups, not as a biological distance.

### 2. Random walk with restart

There were 300 seeds (downstream nodes), with damping 0.85 (restart
probability 0.15) and the original weights.

| Category | Mean score per node | Share of total probability |
|---|---|---|
| Ligands | 0 | 0% |
| Receptors | 4.6–4.9 × 10⁻⁵ | ~5% |
| Downstream nodes | 5.1–5.4 × 10⁻⁵ | ~95% |

Ranges cover different random samples of seeds.

- **No ligand received any probability from any seed.**
- About 95% of the probability stays on downstream nodes (about 15.5% of it
  on the seed itself, because of the restarts) and about 5% reaches
  receptors.
- An even spread over all nodes would give every node 1/19,763 ≈ 5.06 × 10⁻⁵.
  So the share of total probability, not the per-node mean, is what shows
  where the walk goes.

### 3. Reachability analysis

| Weights | Downstream nodes that reach ≥ 1 ligand |
|---|---|
| Original | 0 / 17,674 |
| All = 1 | 0 / 17,674 |

- **No downstream node can reach a ligand.** The 1,005 ligands that are not
  also receptors have no incoming edges in `lr_sig_network`, so no node can
  reach them, by construction of the network.
- Downstream nodes *can* reach receptors: 17,581 of 17,674 downstream nodes
  reach at least one receptor, 1,079 of the 1,084 receptors are reached, and
  each downstream node reaches 1,073 receptors on average.
- If the 282 genes that are both ligand and receptor were kept in the ligand
  group, they would be the only "ligands" reachable.

## Conclusion

The direction of the signal is respected in `lr_sig_network`. Signal started
from downstream nodes can flow to receptors and other downstream nodes, but
never back to ligands. Using all nodes as seed nodes in EV-Net therefore does
not make the signal flow back to ligands.
