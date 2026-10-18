#%%
# =============================================================================
# EXPERIMENT 3 (support) — Toy example of directed reachability in igraph
# -----------------------------------------------------------------------------
# A small mock network to show how Graph.distances(mode="OUT") follows edge
# direction, before applying it to the full NicheNet network in
# experiment3_reachability_to_ligands.py.
#
#   L1 -> R1 -> D1 -> T2        L = ligand, R = receptor,
#         R2 -> T1 <- D2        D = downstream node, T = transcriptional regulator
# =============================================================================
import igraph as ig
import pandas as pd

# Mock edges (from, to, weight)
data = {
    "from": ["L1", "R1", "R2", "D1", "D2"],
    "to":   ["R1", "D1", "T1", "T2", "T1"],
    "weight": [1, 1, 1, 1, 1],   # set all weights = 1
}

mock_network = pd.DataFrame(data)

# Build a directed igraph Graph from the DataFrame
G_mock = ig.Graph.TupleList(
    mock_network.itertuples(index=False, name=None),
    directed=True,
    edge_attrs=["weight"],
)

# Check graph
print(G_mock.summary())
print("Nodes:", G_mock.vs["name"])
print("Edges:", [(e.source, e.target) for e in G_mock.es])

#%%
# Distances following edge direction (mode="OUT")
dist_matrix = G_mock.distances(
    source=["L1", "R1", "D1"],  # sources
    target=["T1", "T2"],        # targets
    weights="weight",
    mode="OUT",
)

print(pd.DataFrame(dist_matrix, index=["L1", "R1", "D1"], columns=["T1", "T2"]))
# Expected:
#      T1   T2
# L1  inf  3.0   L1 -> R1 -> D1 -> T2 (3 steps); no directed path to T1
# R1  inf  2.0
# D1  inf  1.0

# %%
# Going backwards: distance from D1 to the ligand L1
dist = G_mock.distances(
    source=["D1"],
    target=["L1"],
    weights="weight",
    mode="OUT",   # only follow edge direction
)

print("Distance D1 → L1:", dist[0][0])
# Expected: inf. L1 has no incoming edges, so it cannot be reached
# from D1 (or from any other node).

# %%
# =============================================================================
# HOW TO READ THESE METRICS
# -----------------------------------------------------------------------------
# Finite distance (1, 2, 3...): the target can be reached from the source;
#   the value is the number of steps (all weights = 1).
# inf: there is no path that follows the arrows, i.e. the target is NOT
#   reachable. With mode="ALL" (direction ignored) D1 -> L1 would be 2;
#   with mode="OUT" it is inf, which is the behaviour used in experiment 3.
# =============================================================================
