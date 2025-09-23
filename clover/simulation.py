import argparse
import os
import numpy as np
from Bio import Phylo
import sys
import os
import copy
import math
import json
import pandas as pd

# Add parent directory to sys.path
sys.path.append(os.path.abspath(".."))
from simulation_utils import *

def smooth_scores(scores, window=5, poly=2):
    if len(scores) < window:
        return scores
    return savgol_filter(scores, window_length=window, polyorder=poly)

def parsimony_vs_lambda(out, parasites, hosts, host_W_matrices, par_W_matrices, flip_cost_matrix, steps, lower, upper):
    """
    Compute total parsimony for a range of λ values between 0 and 1.
    Returns arrays of lambdas and parsimony scores.
    """
    # filter lambda from 0.5 to 0.9
    lambdas = np.linspace(lower, upper, steps)
    scores = []

    for lam in lambdas:
        cut_result = solve_network_cut(
            out, host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
            flip_cost_matrix=flip_cost_matrix,
            lambda_param=lam
        )
        mat = cut_result["new_matrix"]

        # total parsimony across parasite×host
        total = 0
        for p in parasites:
            leaf_states = {h: cut_result["new_cell_state"][(p,h)] for h in hosts}
            total += sankoff(out["host_trees"][0], leaf_states, np.array([[0,1],[1,0]]))
        for h in hosts:
            leaf_states = {p: cut_result["new_cell_state"][(p,h)] for p in parasites}
            total += sankoff(out["par_trees"][0], leaf_states, np.array([[0,1],[1,0]]))

        scores.append(total)

    return lambdas, scores

def find_elbow_parsimony(out, parasites, hosts, host_W_matrices, par_W_matrices, flip_cost_matrix, lower, upper, outdir="experiments"):
    # Run your sweep
    lambdas, scores = parsimony_vs_lambda(out, parasites, hosts, host_W_matrices, par_W_matrices, flip_cost_matrix, steps=50, lower=lower, upper=upper)

    # Use KneeLocator to detect elbow
    # Smooth scores to reduce noise
    # scores = smooth_scores(scores, window=5, poly=1)
        # Smooth to reduce zig-zag

    lambdas = np.array(lambdas)
    scores = np.array(scores)

    # Smooth to reduce zig-zag
    smooth = smooth_scores(scores, window=3, poly=2)
    smooth_lambdas = lambdas[:len(smooth)]

    # Compute slope
    diffs = np.diff(smooth) / np.diff(smooth_lambdas)

    # Step 1: find steepest slope
    idx_steep = np.argmin(diffs)  # most negative slope
    cutoff_lambda = smooth_lambdas[idx_steep]

    # Step 2: run KneeLocator only on the part after cutoff
    mask_after = lambdas >= cutoff_lambda
    lambdas_after = lambdas[mask_after]
    scores_after = scores[mask_after]

    if len(lambdas_after) < 3:  # fallback if too short
        elbow_lambda = cutoff_lambda
    else:
        kneedle = KneeLocator(
            lambdas_after, scores_after,
            curve="convex", direction="decreasing"
        )
        elbow_lambda = kneedle.elbow if kneedle.elbow else cutoff_lambda

    print(f"[Elbow detection] λ ≈ {elbow_lambda:.3f} (cutoff={cutoff_lambda:.3f})")

    # Plot
    plt.figure(figsize=(7,5))
    plt.plot(lambdas, scores, "bo-", label="Parsimony score")
    plt.axvline(cutoff_lambda, color="orange", linestyle="--",
                label=f"Steep cutoff λ={cutoff_lambda:.3f}")
    if elbow_lambda is not None:
        plt.axvline(elbow_lambda, color="red", linestyle="--",
                    label=f"Elbow λ={elbow_lambda:.3f}")
    plt.xlabel("λ")
    plt.ylabel("Parsimony")
    plt.title("Parsimony vs λ (Steepness + Knee)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "parsimony_vs_lambda_elbow.png"))
    plt.close()

    return elbow_lambda, lambdas, scores


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






def matrix_builder(out, r01_h, r10_h, r01_p, r10_p):
    # Build weight matrices for each parasite’s host tree
    host_W_matrices = []
    for t in out["host_trees"]:
        W, nodes = build_weight_matrix(t, r01_h, r10_h)
        host_W_matrices.append((W, nodes))

    # Build weight matrices for each host’s parasite tree
    par_W_matrices = []
    for t in out["par_trees"]:
        W, nodes = build_weight_matrix(t, r01_p, r10_p)
        par_W_matrices.append((W, nodes))

    # Flip cost matrix
    flip_cost_matrix = build_flip_cost_matrix(len(out["par_leaves"]), len(out["host_leaves"]), cost=1.0)
    return host_W_matrices, par_W_matrices, flip_cost_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--r01_p", type=float, default=0.5)
    parser.add_argument("--r10_p", type=float, default=0.5)
    parser.add_argument("--r01_h", type=float, default=0.5)
    parser.add_argument("--r10_h", type=float, default=0.5)
    parser.add_argument("--host_tree", type=str, required=True, help="Host tree file in Newick format")
    parser.add_argument("--virus_tree", type=str, required=True, help="Virus/parasite tree file in Newick format")
    parser.add_argument("--association_csv", type=str, required=True, help="CSV file with Host,Virus associations")
    parser.add_argument("--corrupt", type=int, default=0)
    parser.add_argument("--lambda_param", type=float, default=0.7, help="Lambda parameter for network cut")
    parser.add_argument("--outdir", type=str, default="experiments", help="Output directory for results")
    parser.add_argument("--original_metrics", type=str, default=None, help="File to save precision/recall/F1 metrics as JSON")
    parser.add_argument("--elbow_metrics", type=str, default=None, help="File to save precision/recall/F1 metrics as JSON")

    args = parser.parse_args()
    r01_p, r10_p, r01_h, r10_h = args.r01_p, args.r10_p, args.r01_h, args.r10_h
    corrupt = args.corrupt
    seed = args.seed
    outdir = args.outdir

    os.makedirs(args.outdir, exist_ok=True)

    # Read associations + trees
    cell_state, base_host_tree, base_par_tree = read_interaction_matrix(
        args.association_csv, args.host_tree, args.virus_tree
    )

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

    # Build the matrix
    mat, parasites, hosts = get_interaction_matrix(out)

    # Save the matrix visualization
    plot_matrix(mat, parasites, hosts, filename=os.path.join(args.outdir, "original_matrix.png"))

    # From here on, same as before...

    host_W_matrices, par_W_matrices, flip_cost_matrix = matrix_builder(out, r01_h, r10_h, r01_p, r10_p)
    
    
    # ---------------------------
    # Step 1: Create corrupted matrix (random 1→0 flips)
    # ---------------------------


    if corrupt > 0:
        random_flips = min(corrupt, np.sum(mat) - 1)
        corrupt_mat, hidden = corrupt_matrix(mat, num_flips=random_flips, seed=seed)

        # Save hidden cells for metrics
        hidden_cells = [(parasites[i], hosts[j]) for i, j in hidden]

        highlight_corrupted = {"corrupted": hidden_cells}
        plot_matrix(corrupt_mat, parasites, hosts,
                    filename=os.path.join(outdir, "corrupted.png"),
                    highlight=highlight_corrupted)

        # Update out["cell_state"] with corrupted values
        corrupt_cell_state = {}
        for i, p in enumerate(parasites):
            for j, h in enumerate(hosts):
                corrupt_cell_state[(p, h)] = int(corrupt_mat[i, j])
        out["cell_state"] = corrupt_cell_state
        

        lambda_param, cut_result = binary_search_lambda(
            out, parasites, hosts, hidden_cells=hidden_cells,
            target_flips=len(hidden_cells),
            host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
            flip_cost_matrix=flip_cost_matrix,
            tol=0, max_iter=20
        )
    else:
        hidden_cells = []
        cut_result = solve_network_cut(
            out, host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
            flip_cost_matrix=flip_cost_matrix,
            lambda_param=args.lambda_param)
        
    # ---------------------------
    # Step 2: Run network cut recovery on corrupted input
    # ---------------------------
    

    new_mat = cut_result["new_matrix"]
    new_cell_state = cut_result["new_cell_state"]
    flips = cut_result["flips"]

    highlight_flipped = {"flipped": [(p, h) for p, h, old, new in flips]}
    print(f"Number of flips performed: {len(flips)}")
    for p, h, old, new in flips:
        print(f"Cell ({p},{h}): {old} -> {new}")

    # ---------------------------
    # Step 3: Save recovered matrix with algorithm flips
    # ---------------------------
    plot_matrix(new_mat, parasites, hosts,
                filename=os.path.join(outdir, "flipped_matrix.png"),
                highlight=highlight_flipped)






    # inside main(), after you have flips and hidden_cells
    metrics = compute_metrics(hidden_cells, flips)
    print("Metrics:", metrics)


    elbow_lambda, lambdas, scores = find_elbow_parsimony(
        out, parasites, hosts, host_W_matrices, par_W_matrices, flip_cost_matrix, lower=0, upper=1, outdir=outdir
    )

    cut_result_elbow = solve_network_cut(
        out, host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
        flip_cost_matrix=flip_cost_matrix,
        lambda_param=elbow_lambda
    )
    flips_elbow = cut_result_elbow["flips"]

    highlight_flipped = {"flipped": [(p, h) for p, h, old, new in flips_elbow]}
    print(f"Number of flips performed (elbow): {len(flips_elbow)}")
    # for p, h, old, new in flips_elbow:
    #     print(f"Cell ({p},{h}): {old} -> {new}")
    plot_matrix(cut_result_elbow["new_matrix"], parasites, hosts,
                filename=os.path.join(outdir, "flipped_matrix_elbow.png"),
                highlight=highlight_flipped)

    metrics_elbow = compute_metrics(hidden_cells, flips_elbow)
    print("Elbow Metrics:", metrics_elbow)

    if args.original_metrics is not None:
        with open(args.original_metrics, "w") as f:
            json.dump(metrics, f, indent=2)
    if args.elbow_metrics is not None:
        with open(args.elbow_metrics, "w") as f:
            json.dump(metrics_elbow, f, indent=2)

if __name__ == "__main__":
    main()
