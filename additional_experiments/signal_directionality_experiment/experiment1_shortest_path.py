#%%
# =============================================================================
# EXPERIMENT 1 — Shortest path analysis
# -----------------------------------------------------------------------------
# Question: are transcriptional regulators (TRs) closer to ligands than to
# "downstream" signaling nodes? If signal flow in the ligand-signaling network
# were strongly hierarchical, ligands should be noticeably farther from (or
# closer to) the TRs than downstream nodes are.
#
# Network: NicheNet mouse lr_sig_network (ligand-receptor + signaling edges),
# exported from R to parquet. Columns: from, to, weight.
#
# Run the cells top to bottom (from this folder). Set USE_UNIT_WEIGHTS below
# to switch between the original NicheNet weights and all weights = 1. The
# histogram is saved to figures/.
# =============================================================================
import random

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# False -> original NicheNet edge weights; True -> every edge has weight 1
# (shortest path = number of steps).
USE_UNIT_WEIGHTS = False

# Fixed seed so the random samples (and therefore the averages) are reproducible
random.seed(42)

#%%
# Load the networks
lr_sig_network = pd.read_parquet("lr_sig_network.parquet")  # graph used below
lr_network = pd.read_parquet("lr_network.parquet")          # ligand -> receptor
sig_network = pd.read_parquet("sig_network.parquet")        # signaling edges

#%%
# -----------------------------------------------------------------------------
# Node categories
# -----------------------------------------------------------------------------
# Ligands: every node that appears as a ligand in the ligand-receptor network
ligands = set(lr_network["from"].unique())

# Transcriptional regulators (TRs): every node that receives a signaling edge.
# Note: this is a broad definition (18,729 nodes, i.e. most of the graph).
targets = set(sig_network["to"].unique())

# Downstream nodes: every node that sends a signaling edge, minus ligands.
# This still overlaps with receptors and with the TRs; experiments 2 and 3
# use a stricter definition.
downstream_nodes = set(sig_network["from"].unique()) - ligands

print("Ligands:", len(ligands))                  # 1,287
print("Targets (TRs):", len(targets))            # 18,729
print("Downstream nodes:", len(downstream_nodes))  # 18,405

#%%
# How much do the categories overlap? (the network is highly interconnected)
print("Ligands ∩ TRs:", len(ligands & targets))
print("Ligands ∩ downstream:", len(ligands & downstream_nodes))
print("TRs ∩ downstream:", len(targets & downstream_nodes))

#%%
# Randomly sample 1,000 sources from each group to keep the computation fast.
# sorted() makes the sample reproducible and is required on Python >= 3.11,
# where random.sample() no longer accepts sets.
sampled_ligands = random.sample(sorted(ligands), 1000)
sampled_downstream_nodes = random.sample(sorted(downstream_nodes), 1000)
targets_list = sorted(targets)

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

print(G_ligand.summary())

#%%
# -----------------------------------------------------------------------------
# Shortest paths (Dijkstra), following edge direction only (mode="OUT")
# -----------------------------------------------------------------------------
# Result: matrix with one row per source and one column per TR. Pairs with
# no directed path get distance inf.
ligand_to_target_paths = np.array(G_ligand.distances(
    source=sampled_ligands,
    target=targets_list,
    weights="weight",
    mode="OUT",
))

downstream_to_target_paths = np.array(G_ligand.distances(
    source=sampled_downstream_nodes,
    target=targets_list,
    weights="weight",
    mode="OUT",
))

#%%
# Average over all reachable (finite) source -> TR pairs
ligand_dists = ligand_to_target_paths[np.isfinite(ligand_to_target_paths)]
downstream_dists = downstream_to_target_paths[np.isfinite(downstream_to_target_paths)]

print("Weights:", "all = 1" if USE_UNIT_WEIGHTS else "original")
print("Average shortest path (ligands -> TRs):", ligand_dists.mean())
print("Average shortest path (downstream -> TRs):", downstream_dists.mean())
print("Fraction of pairs that are reachable (ligands):", np.isfinite(ligand_to_target_paths).mean())
print("Fraction of pairs that are reachable (downstream):", np.isfinite(downstream_to_target_paths).mean())

# Results reported (seed 42 gives almost the same values):
#   original weights: ligands -> TRs 0.22 | downstream -> TRs 0.21
#   weights = 1     : ligands -> TRs 2.41 | downstream -> TRs 2.43

#%%
# Histogram of the two distance distributions
plt.figure(figsize=(8, 5))
plt.hist(ligand_dists, bins=30, alpha=0.6, label="Ligand → Transcriptional regulator")
plt.hist(downstream_dists, bins=30, alpha=0.6, label="Downstream → Transcriptional regulator")
plt.xlabel("Shortest path length")
plt.ylabel("Number of source–target pairs")
plt.title("Distribution of shortest path lengths to transcriptional regulators")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
weights_label = "unit_weights" if USE_UNIT_WEIGHTS else "original_weights"
plt.savefig(f"figures/shortest_path_histogram_{weights_label}.png", dpi=150)
plt.show()

# %%
# =============================================================================
# HOW TO READ THESE METRICS
# -----------------------------------------------------------------------------
# Shortest path length (d): the smallest total edge weight over all directed
#   paths from a source to a target. With all weights = 1 this is simply the
#   number of steps (d = 2 means source -> X -> target).
#
# Original weights: igraph treats weights as COSTS (lengths). In NicheNet a
#   higher weight means a MORE reliable interaction, so the weighted shortest
#   path favours weak edges. Weighted values (~0.2) are only useful for
#   comparing the two groups with each other, not as a biological distance.
#
# Average shortest path: the mean of d over every reachable (finite)
#   source-target pair. Unreachable pairs (inf) are left out; ~0.5% of pairs
#   here.
#
# Interpretation: ligands -> TRs (2.41 steps) and downstream -> TRs (2.43
#   steps) are practically identical. Almost any node reaches almost any TR in
#   2-3 steps ("small-world" network), so distance alone does not separate
#   ligands from downstream nodes.
# =============================================================================
