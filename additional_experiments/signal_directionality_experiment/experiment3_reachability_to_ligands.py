#%%
# =============================================================================
# EXPERIMENT 3 — Reachability analysis
# -----------------------------------------------------------------------------
# Question: following edge direction, can a "downstream" signaling node ever
# reach a ligand? If not, a random walk started from a downstream node cannot
# flow back "upstream" to ligands, so using all nodes of lr_sig_network as
# seed nodes (as EV-Net does) does not reverse the direction of the signal.
#
# Network: NicheNet mouse lr_sig_network (ligand-receptor + signaling edges).
# See experiment3_supportgraph.py for a toy example of the same computation.
#
# Set USE_UNIT_WEIGHTS below to switch between original weights and weights = 1
# (reachability itself does not depend on the weights; only distances do).
# =============================================================================
import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# False -> original NicheNet edge weights; True -> every edge has weight 1
USE_UNIT_WEIGHTS = False

#%%
# Load the networks
lr_sig_network = pd.read_parquet("lr_sig_network.parquet")  # graph used below
lr_network = pd.read_parquet("lr_network.parquet")          # ligand -> receptor
sig_network = pd.read_parquet("sig_network.parquet")        # signaling edges

#%%
# -----------------------------------------------------------------------------
# Node categories
# -----------------------------------------------------------------------------
ligands = set(lr_network["from"].unique())    # 1,287
receptors = set(lr_network["to"].unique())    # 1,084

# 282 nodes act both as ligand and as receptor. They are kept as receptors
# and removed from the ligands, so "ligands" are pure ligands only.
print("Ligands ∩ receptors:", len(ligands & receptors))
ligands = ligands - receptors                 # 1,005

# Downstream nodes: nodes of the signaling network that do not appear
# anywhere in the ligand-receptor network (neither ligands nor receptors).
sig_nodes = set(sig_network["from"]) | set(sig_network["to"])
lr_nodes = set(lr_network["from"]) | set(lr_network["to"])
downstream_nodes = sorted(sig_nodes - lr_nodes)

print("Ligands:", len(ligands))                   # 1,005
print("Receptors:", len(receptors))               # 1,084
print("Downstream nodes:", len(downstream_nodes))  # 17,674
print("Ligands ∩ downstream:", len(ligands & set(downstream_nodes)))  # 0

#%%
# -----------------------------------------------------------------------------
# Build the directed, weighted graph
# -----------------------------------------------------------------------------
G_ligand = ig.Graph.TupleList(
    lr_sig_network.itertuples(index=False, name=None),
    directed=True,
    edge_attrs=["weight"],
)

if USE_UNIT_WEIGHTS:
    G_ligand.es["weight"] = [1.0] * G_ligand.ecount()

#%%
# Sanity check: how many pure ligands have at least one incoming edge?
# If 0, no node can reach them, whatever path is taken.
ligands_in_graph = sorted(ligands & set(G_ligand.vs["name"]))
indegree = dict(zip(G_ligand.vs["name"], G_ligand.indegree()))
print("Pure ligands with incoming edges:", sum(indegree[l] > 0 for l in ligands_in_graph))

#%%
# -----------------------------------------------------------------------------
# Reachability: distances from every downstream node to every ligand
# -----------------------------------------------------------------------------
# mode="OUT" follows edges in their direction only. Distance inf = no
# directed path exists, i.e. the ligand is NOT reachable.
dist_matrix = np.array(G_ligand.distances(
    source=downstream_nodes,
    target=ligands_in_graph,
    weights="weight",
    mode="OUT",
))

#%%
reachable = np.isfinite(dist_matrix)

# For each downstream node: can it reach ANY ligand?
downstream_reachability = reachable.any(axis=1)
reachable_count = downstream_reachability.sum()
total_downstream = len(downstream_nodes)

print("Weights:", "all = 1" if USE_UNIT_WEIGHTS else "original")
print(f"Downstream nodes that reach at least one ligand: {reachable_count} / {total_downstream}")

# Result reported: 0 / 17,674 with original weights and with weights = 1

#%%
# How many ligands each downstream node can reach (all zeros expected)
ligands_reached_per_node = reachable.sum(axis=1)
print("Average ligands reachable per downstream node:", ligands_reached_per_node.mean())

plt.hist(ligands_reached_per_node, bins=30)
plt.xlabel("Number of ligands reachable")
plt.ylabel("Number of downstream nodes")
plt.show()

#%%
# -----------------------------------------------------------------------------
# For comparison: reachability of receptors
# -----------------------------------------------------------------------------
# Unlike ligands, receptors do have incoming signaling edges, so downstream
# nodes can reach them.
receptors_in_graph = sorted(receptors & set(G_ligand.vs["name"]))
reachable_rec = np.isfinite(np.array(G_ligand.distances(
    source=downstream_nodes,
    target=receptors_in_graph,
    weights="weight",
    mode="OUT",
)))
print(f"Downstream nodes that reach at least one receptor: {reachable_rec.any(axis=1).sum()} / {total_downstream}")
print("Average receptors reachable per downstream node:", reachable_rec.sum(axis=1).mean())
print("Distinct receptors reached:", reachable_rec.any(axis=0).sum(), "/", len(receptors_in_graph))

# %%
# =============================================================================
# HOW TO READ THESE METRICS
# -----------------------------------------------------------------------------
# Reachable: node B is reachable from node A if there is at least one path
#   A -> ... -> B that follows edge directions. igraph returns a finite
#   distance if B is reachable and inf if not. Weights only change the
#   distance value, never whether B is reachable, which is why original
#   weights and weights = 1 give the same counts.
#
# "Downstream nodes that reach at least one ligand": the number of downstream
#   nodes with at least one finite distance to a ligand. The result is 0.
#
# Why it is 0: pure ligands have NO incoming edges in lr_sig_network (see
#   the sanity check above). The network is built so that ligands are only
#   sources, so no walk can ever arrive at a ligand. This holds by
#   construction for every node, not only for the downstream nodes.
#
# Receptors: downstream nodes DO reach receptors (almost all of the 1,084),
#   because receptors also receive signaling edges. The 282 nodes that are
#   both ligand and receptor are reachable for the same reason; they are the
#   only "ligands" reachable when receptors are not removed from the ligands.
# =============================================================================
