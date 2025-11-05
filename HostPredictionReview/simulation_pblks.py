import pandas as pd
import os
import argparse
import numpy as np
from Bio import Phylo
import sys
import copy
import math
import json
import matplotlib.pyplot as plt

# Add parent directory to sys.path
sys.path.append(os.path.abspath(".."))
import importlib
# importlib.reload(sys.modules['simulation_utils'])
from simulation_utils import *


def read_interaction_matrix(csv_file, host_tree_file, virus_tree_file):
    """
    Read host-virus interaction matrix CSV and reorder according to tree leaves.
    Treat viruses as parasites.
    
    Returns:
        cell_state: dict with keys (virus, host) -> 0/1
        host_tree: Bio.Phylo tree object
        virus_tree: Bio.Phylo tree object
    """

    # Load trees
    host_tree = Phylo.read(host_tree_file, "newick")
    virus_tree = Phylo.read(virus_tree_file, "newick")

    # Get leaf order
    host_order = [t.name for t in host_tree.get_terminals()]
    virus_order = [t.name for t in virus_tree.get_terminals()]

    print(f"Number of hosts in tree: {len(host_order)}")
    print(f"Number of viruses in tree: {len(virus_order)}")

    # Load interaction matrix CSV
    df = pd.read_csv(csv_file, index_col=0)  # assuming rows=viruses, columns=hosts
    # Ensure names match the expected tree leaf names
    df = df.loc[df.index.intersection(virus_order), df.columns.intersection(host_order)]

    # Reorder according to tree leaves
    df = df.reindex(index=virus_order, columns=host_order, fill_value=0)

    # Build cell_state dictionary
    cell_state = {}
    for virus in virus_order:
        for host in host_order:
            cell_state[(virus, host)] = int(df.at[virus, host])

    return cell_state, host_tree, virus_tree

def parse_phist_costs(csv_path, hosts, parasites):
    """
    Build a flip cost matrix using k-mer CSV file, mapping parasites (rows) × hosts (columns).
    Each cell = k-mer count from CSV, or 0 if missing.
    """
    with open(csv_path) as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Extract phage names from the header (3rd column onwards)
    header_line = next(l for l in lines if l.startswith("kmer-length"))
    header_parts = header_line.split(",")
    phage_names = [p.replace(".fa", "") for p in header_parts[2:] if p.startswith("k141_")]
    
    # Map phage index (1-based) → phage name
    phage_index = {i+1: name for i, name in enumerate(phage_names)}
    
    # Prepare empty cost matrix: rows = parasites, cols = hosts
    cost_matrix = np.zeros((len(parasites), len(hosts)), dtype=int)
    
    for line in lines:
        if not line.startswith("bin."):
            continue
        
        parts = line.split(",")
        host_name = parts[0].replace(".fa", "")
        if host_name not in hosts:
            continue
        
        host_idx = hosts.index(host_name)
        
        # Iterate over all kmer pairs like 4:10
        for token in parts[2:]:
            if ":" not in token:
                continue
            try:
                phage_num, value = map(int, token.split(":"))
            except ValueError:
                continue
            
            if phage_num not in phage_index:
                continue
            phage_name = phage_index[phage_num]
            
            if phage_name not in parasites:
                continue
            
            phage_idx = parasites.index(phage_name)
            cost_matrix[phage_idx, host_idx] = value

    # row_sums = cost_matrix.sum(axis=1, keepdims=True)
    #     # Avoid division by zero
    # row_sums[row_sums == 0] = 1
    # cost_matrix = cost_matrix / row_sums
            

    # Return DataFrame (rows = viruses, columns = hosts)

    return cost_matrix

def parse_pblks_costs(csv_path, hosts, parasites):
    """
    Parse a binary virus–host matrix and assign weighted costs.
    For each row (virus):
      - If entry == 1: cost = column_index / (number of 1s in that row)
      - If entry == 0: cost = 1
    Column priority = left to right (1-based index).
    """
    # Load CSV
    df = pd.read_csv(csv_path, index_col=0)

    # Clean up names
    df.index = df.index.str.replace(".fa", "", regex=False)
    df.columns = df.columns.str.replace(".fa", "", regex=False)

    # Keep only intersecting rows/cols before weighting
    df = df.loc[df.index.intersection(parasites), df.columns.intersection(hosts)]

    # Convert to numeric (in case CSV parsed as str)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

    # Initialize cost matrix (float)
    cost_df = pd.DataFrame(2.0, index=df.index, columns=df.columns)

    # Compute weighted costs
    for i, row in df.iterrows():
        ones = row.values == 1
        num_ones = ones.sum()
        if num_ones == 0:
            continue  # all costs remain 1
        col_indices = np.arange(1, len(row) + 1)  # 1-based column indices
        weights = col_indices / len(df.columns)
        cost_df.loc[i, ones] = weights[ones]

    # Reindex to match provided order and fill missing with 1
    cost_df = cost_df.reindex(index=parasites, columns=hosts, fill_value=2.0)

    # Return as numpy float array
    return cost_df.to_numpy(dtype=float)


def parse_wish_costs(matrix_path, hosts, parasites):
    """
    Parse a virus–host cost matrix and assign rank-based costs
    derived from log-likelihood scores.

    For each phage (column):
      - Rank hosts by descending score (higher = better)
      - Assign cost = rank / total_hosts (so smaller = better)

    Parameters
    ----------
    matrix_path : str
        Path to the .matrix file (TSV: rows = hosts, columns = parasites)
    hosts : list[str]
        List of host (MAG/bin) names to retain and order.
    parasites : list[str]
        List of parasite (phage) names to retain and order.

    Returns
    -------
    cost_matrix : np.ndarray
        2D array (len(parasites) × len(hosts)), ordered as given.
    """

    # Load tab-separated matrix (first column = row names)
    df = pd.read_csv(matrix_path, sep='\t', index_col=0)

    # Clean up names
    df.index = df.index.str.replace(".fa", "", regex=False)
    df.columns = df.columns.str.replace(".fa", "", regex=False)

    # Filter to intersection first
    df = df.loc[df.index.intersection(hosts), df.columns.intersection(parasites)]

    # Reindex to given order (fill missing with NaN)
    df = df.reindex(index=hosts, columns=parasites).astype(float)

    # Initialize cost DataFrame
    cost_df = pd.DataFrame(index=hosts, columns=parasites, dtype=float)

    # Assign rank-based costs for each parasite (phage)
    for phage in df.columns:
        scores = df[phage]
        n_hosts = scores.notna().sum()
        if n_hosts == 0:
            cost_df[phage] = np.nan
            continue

        # Rank hosts: high score → low rank (1 = best)
        ranked_hosts = scores.rank(ascending=False, method='first')-1

        # Cost = rank / total_hosts (so 1st = 1/n, last = 1)
        cost_df[phage] = ranked_hosts / n_hosts

    # Transpose → rows = parasites, columns = hosts
    cost_df = cost_df.T.reindex(index=parasites, columns=hosts, fill_value=1.0)

    # Return as NumPy array
    return cost_df.to_numpy(dtype=float)


def parse_iphop_costs(csv_path, ground_truth_path, hosts, parasites):
    """
    Build a cost matrix from ranked host predictions per phage.

    Each phage's costs are assigned relatively:
      - Predicted hosts (found in host list) ranked 1, 2, 3, ...
      - All other hosts assigned (max_rank + 1).

    Parameters
    ----------
    csv_path : str
        Path to CSV with 'phage' and 'MAGs' (comma-separated list sorted by confidence).
    ground_truth_path : str
        Path to ground truth TSV with 'Phage' and 'MAGs' columns.
    hosts : list[str]
        Ordered list of host names (columns).
    parasites : list[str]
        Ordered list of phage names (rows).

    Returns
    -------
    np.ndarray
        Cost matrix (float) with shape [len(parasites), len(hosts)].
        Lower = better (1 = top prediction, larger = worse).
    """

    # --- Load and parse ---
    pred = pd.read_csv(csv_path)
    truth = pd.read_csv(ground_truth_path, sep="\t")

    # Build lookup for true MAG per phage (not directly used for cost, but could be later)
    true_host_map = {row["Phage"]: row["MAGs"] for _, row in truth.iterrows()}

    # Initialize cost DataFrame
    cost_df = pd.DataFrame(np.nan, index=parasites, columns=hosts, dtype=float)

    # Iterate through each phage in predictions
    for _, row in pred.iterrows():
        phage = row["phage"]
        if phage not in cost_df.index:
            continue

        # Parse predicted hosts
        pred_hosts = [h.strip() for h in str(row["MAGs"]).split(",") if h.strip()]

        # Keep only those that exist in the host list
        matched_hosts = [h for h in pred_hosts if h in hosts]

        if len(matched_hosts) == 0:
            # No valid hosts found — assign all hosts max rank = len(hosts)
            cost_df.loc[phage, :] = float(len(hosts))
            continue

        # Assign relative ranks (1, 2, 3, ...)
        for rank, host in enumerate(matched_hosts, start=1):
            cost_df.loc[phage, host] = float(rank)

        # For all other hosts, assign (max_rank + 1)
        max_rank = len(matched_hosts)
        cost_df.loc[phage, cost_df.columns.difference(matched_hosts)] = float(max_rank + 1)

    # Fill any remaining NaN (for phages missing predictions) with len(hosts) + 1
    cost_df = cost_df.fillna(float(len(hosts)))

    # normalize each row by number of hosts
    for phage in cost_df.index:
        row = cost_df.loc[phage]
        num_hosts = len(row)
        cost_df.loc[phage] = row / num_hosts

    # Reindex for exact order
    cost_df = cost_df.reindex(index=parasites, columns=hosts)

    return cost_df.to_numpy(dtype=float)

def matrix_builder(out, r01_h, r10_h, r01_p, r10_p):
    host_W_matrices = []
    for t in out["host_trees"]:
        W, nodes = build_weight_matrix(t, r01_h, r10_h, scale=100.0)
        # scale W
        # W[W!=0] =  np.inf
        W = W*100
        host_W_matrices.append((W, nodes))

    # Build weight matrices for each host’s parasite tree
    par_W_matrices = []
    for t in out["par_trees"]:
        W, nodes = build_weight_matrix(t, r01_p, r10_p, scale=100.0)
        # scale W
        W = W*100
        par_W_matrices.append((W, nodes))
    # Flip cost matrix
    # flip_cost_matrix = build_flip_cost_matrix(len(out["par_leaves"]), len(out["host_leaves"]), cost=100.0)
    # flip_cost_matrix=-np.log(flipping_cost_matrix)
    return host_W_matrices, par_W_matrices

def flip_specific_cells(mat, parasites, hosts, flip_targets):
    """
    Flip given cells (1->0) based on Phage-HOST pairs.
    Also set flip cost rules:
      - Flipped cells -> cost = 0.5
      - All other phage rows (not in flip targets) -> inf for all hosts
    
    Args:
        mat: numpy array, original matrix
        parasites: list of row labels (Phages)
        hosts: list of column labels (MAGs)
        flip_targets: list of tuples (Phage, MAGs) to flip
        flip_cost_matrix: numpy array (same shape as mat) to modify
    
    Returns:
        mat_corrupt: numpy array with specified cells flipped
        flipped_indices: list of (i,j) indices that were flipped
        flip_cost_matrix: modified cost matrix
    """
    mat_corrupt = mat.copy()
    flipped_indices = []

    parasite_idx = {p: i for i, p in enumerate(parasites)}
    host_idx = {h: j for j, h in enumerate(hosts)}

    # Track which phages are in flip targets
    targeted_phages = set()

    for phage, host in flip_targets:
        if phage in parasite_idx and host in host_idx:
            i, j = parasite_idx[phage], host_idx[host]
            # flip_cost_matrix[i, :] = 50
            targeted_phages.add(phage)
            if mat_corrupt[i, j] == 1:
                mat_corrupt[i, j] = 0
                # flip_cost_matrix[i, j] = 0.0
                flipped_indices.append((i, j))

    # Set inf for all non-targeted phage rows
    for p in parasites:
        if p not in targeted_phages:
            i = parasite_idx[p]
            # flip_cost_matrix[i, :] = 1
            mat_corrupt[i, :] = 1

    return mat_corrupt, flipped_indices

# --- datasets ---
datasets = [
    {
        "name": "cow",
        "host_tree": "host_taxonomy_tree_cow.newick",
        "virus_tree": "virus_tree_cow.newick",
        "interaction_tsv": "cow_dataset_hic_gtdb.tsv"
    },
    {
        "name": "gut",
        "host_tree": "host_taxonomy_tree_gut.newick",
        "virus_tree": "virus_tree_human_gut.newick",
        "interaction_tsv": "gut_dataset_hic_gtdb.tsv"
    },
    {
        "name": "water",
        "host_tree": "host_taxonomy_tree_wastewater.newick",
        "virus_tree": "virus_tree_wastewater.newick",
        "interaction_tsv": "water_dataset_hic_gtdb.tsv"
    }
]

# --- target flips to test ---
target_flips_list = [1, 2, 3, 4, 5, 10]

# --- output dir ---
os.makedirs("experiments", exist_ok=True)


def run_experiment(dataset_name, host_tree_file, virus_tree_file, interaction_tsv, target_flips):
    """
    Runs your full pipeline for a given dataset and number of target flips.
    Returns (accuracy, numerator, denominator)
    """
    print(f"\n🔬 Running dataset={dataset_name}, target_flips={target_flips}")

    # Build interaction matrix
    output_csv = f"interaction_matrix_{dataset_name}.csv"
    df = pd.read_csv(interaction_tsv, sep="\t")

    phages = df['Phage'].unique()
    hosts = df['MAGs'].unique()

    host_tree = Phylo.read(host_tree_file, "newick")
    virus_tree = Phylo.read(virus_tree_file, "newick")

    ordered_hosts = [t.name for t in host_tree.get_terminals() if t.name in hosts]
    ordered_phages = [t.name for t in virus_tree.get_terminals() if t.name in phages]

    interaction_matrix = pd.DataFrame(0, index=ordered_phages, columns=ordered_hosts)
    for _, row in df.iterrows():
        if row['Phage'] in interaction_matrix.index and row['MAGs'] in interaction_matrix.columns:
            interaction_matrix.loc[row['Phage'], row['MAGs']] = 1

    interaction_matrix.to_csv(output_csv)

    # --- Read and build experiment state ---
    cell_state, base_host_tree, base_par_tree = read_interaction_matrix(output_csv, host_tree_file, virus_tree_file)
    base_host_tree = rescale_tree(base_host_tree)
    base_par_tree = rescale_tree(base_par_tree)

    host_leaves = [t.name for t in base_host_tree.get_terminals()]
    par_leaves = [t.name for t in base_par_tree.get_terminals()]

    out = {
        "host_trees": [copy.deepcopy(base_host_tree) for _ in range(len(par_leaves))],
        "par_trees": [copy.deepcopy(base_par_tree) for _ in range(len(host_leaves))],
        "cell_state": cell_state,
        "host_leaves": host_leaves,
        "par_leaves": par_leaves
    }

    mat, parasites, hosts = get_interaction_matrix(out)

    # --- cost matrix (example: WIsH, modify if you want PHIST etc.) ---
    flip_cost_matrix = parse_pblks_costs(f"../../HostPredictionReview/Benchmark/Task2-MetaHiC/Results/pblks/{dataset_name}.csv", hosts, parasites)
    flip_cost_matrix = flip_cost_matrix*100
    # --- weights ---
    r01_p, r10_p, r01_h, r10_h = 0.7, 0.3, 0.5, 0.5
    host_W_matrices, par_W_matrices = matrix_builder(out, r01_h, r10_h, r01_p, r10_p)

    # --- run for each phage–host pair ---
    total = 0
    acc = 0
    for _, row in df.iterrows():
        phage, mag = row["Phage"], row["MAGs"]
        flip_targets = [(phage, mag)]

        corrupt_mat, hidden = flip_specific_cells(mat, parasites, hosts, flip_targets)
        hidden_cells = [(parasites[i], hosts[j]) for i, j in hidden]

        corrupt_cell_state = {(p, h): int(corrupt_mat[i, j])
                              for i, p in enumerate(parasites)
                              for j, h in enumerate(hosts)}
        out["cell_state"] = corrupt_cell_state

        lambda_param, cut_result = binary_search_lambda(
            out,
            parasites,
            hosts,
            hidden_cells=hidden_cells,
            target_flips=target_flips,
            host_W_matrices=host_W_matrices,
            par_W_matrices=par_W_matrices,
            flip_cost_matrix=flip_cost_matrix,
            tol=0,
            max_iter=20,
        )

        flips = cut_result["flips"]
        metrics = compute_metrics(hidden_cells, flips)

        total += 1
        if metrics["recall"] == 1:
            acc += 1

    accuracy = acc / total if total > 0 else 0
    print(f"✅ {dataset_name}, flips={target_flips} → Accuracy {acc}/{total} = {accuracy:.2f}")
    return accuracy, acc, total


# --- MAIN LOOP ---
summary_rows = []
for ds in datasets:
    row = {"dataset": ds["name"]}
    for n in target_flips_list:
        acc, num, den = run_experiment(ds["name"], ds["host_tree"], ds["virus_tree"], ds["interaction_tsv"], n)
        row[f"{n}_hit"] = f"{acc:.2f}"
        row[f"{n}_hit_numden"] = f"{num}/{den}"
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("experiments/flip_accuracy_summary_pblks.tsv", sep="\t", index=False)
print("\n📊 Summary saved to experiments/flip_accuracy_summary_pblks.tsv")
print(summary_df)
