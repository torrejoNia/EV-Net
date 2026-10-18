#%%
# =============================================================================
# EXPERIMENT 2 — Random walk with restart (RWR)
# -----------------------------------------------------------------------------
# Question: if a random walk with restart (personalized PageRank) starts from
# a "downstream" signaling node, does the probability flow back upstream to
# ligands/receptors, or does it stay downstream (toward the regulators)?
# This is the situation EV-Net creates when all nodes of lr_sig_network are
# used as seed nodes.
#
# Network: NicheNet mouse lr_sig_network (ligand-receptor + signaling edges),
# original weights.
# =============================================================================
import random

import igraph as ig
import numpy as np
import pandas as pd

# Fixed seed so the random samples are reproducible
random.seed(42)

#%%
# Load the networks
lr_sig_network = pd.read_parquet("lr_sig_network.parquet")  # graph used below
lr_network = pd.read_parquet("lr_network.parquet")          # ligand -> receptor
sig_network = pd.read_parquet("sig_network.parquet")        # signaling edges

#%%
# -----------------------------------------------------------------------------
# Node categories (same definitions as experiment 3)
# -----------------------------------------------------------------------------
ligands = set(lr_network["from"].unique())
receptors = set(lr_network["to"].unique())

# Downstream nodes: nodes of the signaling network that do not appear
# anywhere in the ligand-receptor network (neither ligands nor receptors).
sig_nodes = set(sig_network["from"]) | set(sig_network["to"])
lr_nodes = set(lr_network["from"]) | set(lr_network["to"])
downstream_nodes = sorted(sig_nodes - lr_nodes)

#%%
# Overlaps between categories
print("Ligands ∩ receptors:", len(ligands & receptors))                  # 282
print("Ligands ∩ downstream:", len(ligands & set(downstream_nodes)))      # 0
print("Receptors ∩ downstream:", len(receptors & set(downstream_nodes)))  # 0

#%%
# Nodes that are both ligand and receptor are counted as receptors only
ligands = ligands - receptors

print(f"Ligands: {len(ligands)}")                    # 1,005
print(f"Receptors: {len(receptors)}")                # 1,084
print(f"Downstream nodes: {len(downstream_nodes)}")  # 17,674

# %%
# Two independent random samples of downstream nodes:
# - sampled_downstream_nodes (1,050): the "downstream" category whose scores
#   are averaged below
# - query_downstream_nodes (300): the seeds, one RWR per seed
sampled_downstream_nodes = set(random.sample(downstream_nodes, 1050))
query_downstream_nodes = random.sample(downstream_nodes, 300)

#%%
# -----------------------------------------------------------------------------
# Build the directed, weighted graph
# -----------------------------------------------------------------------------
G = ig.Graph.TupleList(
    lr_sig_network.itertuples(index=False, name=None),
    directed=True,
    edge_attrs=["weight"],
)
print(G.summary())

#%%
def run_ppr(seed_nodes, G, damping=0.85, cutoff=None):
    """Random walk with restart (personalized PageRank) from seed_nodes.

    At every step the walker follows an outgoing edge (chosen proportionally
    to its weight) with probability `damping`, or jumps back to a seed node
    with probability 1 - damping. Returns one score per vertex of G; the
    scores sum to 1.
    """
    # Convert seed node names -> vertex indices
    seed_indices = [G.vs.find(name=s).index for s in seed_nodes]

    # Restart (reset) vector: uniform probability over the seed nodes
    reset = np.zeros(G.vcount())
    for idx in seed_indices:
        reset[idx] = 1.0 / len(seed_indices)

    scores = np.array(G.personalized_pagerank(
        reset=reset.tolist(),  # must be a list, not a dict
        damping=damping,
        directed=True,         # follow edge direction
        weights="weight",
    ))

    # Optional: set scores below `cutoff` to 0 (not used for the results)
    if cutoff is not None:
        scores[scores < cutoff] = 0

    return scores

#%%
# Category sets used to summarize the scores
ligands_set = set(ligands)
receptors_set = set(receptors)
downstream_set = set(sampled_downstream_nodes)

#%%
# -----------------------------------------------------------------------------
# One RWR per seed (300 runs): rows = seeds, columns = all graph nodes
# -----------------------------------------------------------------------------
all_scores = []

for dn in query_downstream_nodes:
    scores = run_ppr([dn], G, damping=0.85, cutoff=None)
    all_scores.append(scores)

scores_df = pd.DataFrame(all_scores, columns=G.vs["name"])

#%%
# Average RWR score per node of each category (mean over seeds and nodes)
is_lig = scores_df.columns.isin(ligands_set)
is_rec = scores_df.columns.isin(receptors_set)
is_down = scores_df.columns.isin(downstream_set)

category_means = {
    "ligand": scores_df.loc[:, is_lig].mean().mean(),
    "receptor": scores_df.loc[:, is_rec].mean().mean(),
    "downstream": scores_df.loc[:, is_down].mean().mean(),
    "other": scores_df.loc[:, ~(is_lig | is_rec | is_down)].mean().mean(),
}
print(category_means)
print("Uniform baseline (1 / number of nodes):", 1 / G.vcount())

# Result reported (original weights):
#   ligand 0 | receptor 4.59e-05 | downstream 5.23e-05
#   (uniform baseline 1/19,763 = 5.06e-05)

#%%
# Fraction of the total probability that ends up in each category (per seed,
# then averaged). This shows WHERE the walk goes, which the per-node means
# above cannot, because the categories differ a lot in size.
is_all_down = scores_df.columns.isin(set(downstream_nodes))
mass = {
    "ligand": scores_df.loc[:, is_lig].sum(axis=1).mean(),
    "receptor": scores_df.loc[:, is_rec].sum(axis=1).mean(),
    "all downstream nodes (17,674)": scores_df.loc[:, is_all_down].sum(axis=1).mean(),
}
print(mass)
print("Max score of any ligand over all seeds:", scores_df.loc[:, is_lig].to_numpy().max())

# Re-run with seed 42: ligand 0.00 | receptor ~0.05 | downstream ~0.95

# %%
print("Ligands in graph:", len(ligands_set & set(scores_df.columns)))
print("Receptors in graph:", len(receptors_set & set(scores_df.columns)))
print("Downstream in graph:", len(downstream_set & set(scores_df.columns)))

# %%
# =============================================================================
# HOW TO READ THESE METRICS
# -----------------------------------------------------------------------------
# RWR score of a node: the long-run probability of finding the walker on that
#   node. The walker starts at the seed, follows outgoing edges (proportionally
#   to their weights) with probability 0.85, and jumps back to the seed with
#   probability 0.15. The scores of all 19,763 nodes add up to 1 for each seed.
#
# Mean score per category: the average score per node of the category, over
#   all 300 seeds. Because the scores sum to 1, a perfectly even spread would
#   give every node 1/19,763 = 5.06e-05. Receptor (4.59e-05) and downstream
#   (5.23e-05) means are close to this value, so these two numbers alone do
#   not say that the walk "stays" in those groups.
#
# Ligand score = 0: the important result. Not a single ligand gets any
#   probability from any of the 300 seeds. Ligands have no incoming edges in
#   lr_sig_network (see experiment 3), so a walker can never step onto one.
#   The only way to reach a ligand is to restart there, which happens only if
#   the ligand is itself a seed.
#
# Fraction of probability per category: the sum of the scores of all nodes
#   in a category. This is the right measure for "where does the diffusion
#   go". Re-run with seed 42: ligands 0%, receptors ~5%, downstream nodes
#   ~95% (of which ~15.5% stays on the seed itself, from the restarts).
# =============================================================================
