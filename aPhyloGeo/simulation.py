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


def read_association_csv(csv_file, host_tree_file, parasite_tree_file):
    """Read host-virus associations and build cell_state + hidden_cells."""

    # Load trees
    host_tree = Phylo.read(host_tree_file, "newick")
    parasite_tree = Phylo.read(parasite_tree_file, "newick")

    # Host tree figure
    fig, ax = plt.subplots(figsize=(20, 15))
    Phylo.draw(host_tree, do_show=False, axes=ax)
    host_tree_png = host_tree_file.replace(".nwk", ".png").replace(".newick", ".png")
    plt.savefig(host_tree_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Parasite tree figure
    fig, ax = plt.subplots(figsize=(20, 15))
    Phylo.draw(parasite_tree, do_show=False, axes=ax)
    parasite_tree_png = parasite_tree_file.replace(".nwk", ".png").replace(".newick", ".png")
    plt.savefig(parasite_tree_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    host_order = [t.name for t in host_tree.get_terminals()]
    parasite_order = [t.name for t in parasite_tree.get_terminals()]

    print(f"Number of hosts in tree: {len(host_order)}")
    print(f"Number of parasites in tree: {len(parasite_order)}")

    # Load associations
    df = pd.read_csv(csv_file)
    associations_raw = [(row["Virus"].strip(), row["Host"].strip()) for _, row in df.iterrows()]

    associations = set()
    duplicates = []
    skipped = []

    for virus, host in associations_raw:
        if virus not in parasite_order or host not in host_order:
            skipped.append((host, virus))
            continue
        if (virus, host) in associations:
            duplicates.append((host, virus))
        associations.add((virus, host))

    # Build cell_state: 1 if association exists, 0 otherwise
    cell_state = {}
    hidden_cells = []  # nothing hidden here, keep empty
    

    for p in parasite_order:
        for h in host_order:
            cell_state[(p, h)] = 1 if (p, h) in associations else 0

    # # count number of 1s
    # num_ones = sum(cell_state.values())

    # # ---- Debug Info ----
    # if skipped:
    #     print("\n⚠️ Skipped associations (not in trees):")
    #     for h, v in skipped:
    #         print(f"  Host={h}, Virus={v}")

    # if duplicates:
    #     print("\nℹ️ Duplicate associations found (already counted once):")
    #     for h, v in duplicates:
    #         print(f"  Host={h}, Virus={v}")

    # print(f"\n✅ Associations kept: {len(associations)}")
    # print(f"❌ Associations skipped: {len(skipped)}")
    # print(f"🔁 Duplicates ignored: {len(duplicates)}")
    # print(f"Total associations (1s) in the matrix: {num_ones}\n")

    return cell_state, host_tree, parasite_tree


def get_top_interactions_from_parafit(parafit_file, association_csv, top_n=30):
    """
    Read ParaFit cailliez results and return top N interactions with host and virus names.
    
    Args:
        parafit_file: Path to output_parafit_cailliez.txt file
        association_csv: Path to association CSV file 
        top_n: Number of top interactions to return (default 30)
    
    Returns:
        List of tuples: [(virus_id, host_full_name), ...]
    """
    
    # Host and virus mappings from the R analysis (1-based indexing)
    host_mapping = {
        1: 'ON640726', 2: 'AB085735', 3: 'ON640722', 4: 'AB085738', 5: 'AB085739', 
        6: 'AB106608', 7: 'AB085731', 8: 'MZ936290', 9: 'KP972690', 10: 'HM134917', 
        11: 'ON012504', 12: 'KX261916', 13: 'JX502551', 14: 'OP894116', 15: 'DQ888677', 
        16: 'ON640662', 17: 'MN366288'
    }
    
    virus_mapping = {
        1: 'EF065505', 2: 'KJ473822', 3: 'EF065509', 4: 'KJ473820', 5: 'KJ473821', 
        6: 'NC', 7: 'GU190215', 8: 'MG772933', 9: 'MG772934', 10: 'MN996532', 
        11: 'DQ022305', 12: 'DQ084200', 13: 'GQ153545', 14: 'GQ153544', 15: 'GQ153546', 
        16: 'GQ153548', 17: 'GQ153540', 18: 'GQ153539', 19: 'GQ153541', 20: 'DQ084199', 
        21: 'GQ153547', 22: 'GQ153543', 23: 'GQ153542', 24: 'DQ648857', 25: 'DQ412043', 
        26: 'KJ473814', 27: 'KJ473812', 28: 'KJ473813', 29: 'KJ473811', 30: 'JX993987', 
        31: 'KY770859', 32: 'KY770858', 33: 'KJ473816', 34: 'KT444582', 35: 'KY417150', 
        36: 'KY417146', 37: 'KY417144', 38: 'KC881005', 39: 'KF367457', 40: 'KC881006', 
        41: 'KY417152', 42: 'KY417151', 43: 'MK211376', 44: 'KY417147', 45: 'KY417148', 
        46: 'KY417143', 47: 'KY417149', 48: 'MK211377', 49: 'KY417142', 50: 'FJ588686', 
        51: 'DQ071615', 52: 'KJ473815', 53: 'KF569996', 54: 'JX993988', 55: 'KU973692', 
        56: 'KF636752', 57: 'KJ473795', 58: 'EU420138', 59: 'KJ473796', 60: 'EU420137', 
        61: 'EU420139', 62: 'KJ473797', 63: 'KJ473800', 64: 'KJ473798', 65: 'KJ473799', 
        66: 'KJ473807', 67: 'KJ473806', 68: 'EF203065', 69: 'KJ473808'
    }
    
    # Read associations to get full host names
    df_assoc = pd.read_csv(association_csv)
    host_id_to_full_name = {}
    for _, row in df_assoc.iterrows():
        host_full = row["Host"].strip()
        # Extract just the accession ID part (before the underscore if present)
        host_id = host_full.split('.')[0]
        host_id_to_full_name[host_id] = host_full
    
    # Read ParaFit results
    interactions = []
    with open(parafit_file, 'r') as f:
        lines = f.readlines()
    
    # Find the start of individual test results
    start_line = None
    for i, line in enumerate(lines):
        if "Host Parasite    F1.stat" in line:
            start_line = i + 1
            break
    
    if start_line is None:
        raise ValueError("Could not find ParaFit results in file")
    
    # Parse results and extract F1 statistics
    for line in lines[start_line:]:
        line = line.strip()
        if not line or line.startswith('Number of'):
            break
            
        # Parse line like: " [1,]    1       58  210.06975 0.001  0.0052777042 0.001"
        parts = line.split()
        if len(parts) >= 4 and parts[0].startswith('['):
            try:
                host_id = int(parts[1])
                virus_id = int(parts[2])
                f1_stat = float(parts[3])
                p_value = float(parts[4])
                
                # Only include significant interactions (p <= 0.05)
                if p_value <= 0.05:
                    interactions.append((host_id, virus_id, f1_stat, p_value))
            except (ValueError, IndexError):
                continue
    
    # Sort by F1 statistic (descending) and take top N
    interactions.sort(key=lambda x: x[2], reverse=True)
    top_interactions = interactions[:top_n]
    
    # Convert to required format: (virus_id, host_full_name)
    result = []
    for host_id, virus_id, f1_stat, p_value in top_interactions:
        # Get virus name
        virus_name = virus_mapping.get(virus_id, f"VIRUS_{virus_id}")
        
        # Get host ID and then full name
        host_short_id = host_mapping.get(host_id, f"HOST_{host_id}")
        host_full_name = host_id_to_full_name.get(host_short_id, host_short_id)
        
        result.append((virus_name, host_full_name))
    return result


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
    cell_state, base_host_tree, base_par_tree = read_association_csv(
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
        top_interactions = get_top_interactions_from_parafit(
            parafit_file="output_parafit_cailliez.txt",
            association_csv=args.association_csv,
            top_n=69
        )
        
        forced_hidden_cells = top_interactions[0:5]
        
        hidden_cells = []
        for (p, h) in forced_hidden_cells:
            if (p, h) in out["cell_state"] and out["cell_state"][(p, h)] == 1:
                out["cell_state"][(p, h)] = 0
                hidden_cells.append((p, h))

        # Build new corrupt_mat from updated cell_state
        corrupt_mat = np.zeros((len(parasites), len(hosts)), dtype=int)
        for i, p in enumerate(parasites):
            for j, h in enumerate(hosts):
                corrupt_mat[i, j] = out["cell_state"][(p, h)]

        # Plot corrupted matrix with highlights
        highlight_corrupted = {"corrupted": hidden_cells}
        plot_matrix(corrupt_mat, parasites, hosts,
                    filename=os.path.join(outdir, "corrupted.png"),
                    highlight=highlight_corrupted)
        
        
        if len(hidden_cells) == 0:
            cut_result = solve_network_cut(
                out, host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
                flip_cost_matrix=flip_cost_matrix,
                lambda_param=args.lambda_param)
        
        else:

            lambda_param, cut_result = binary_search_lambda(
                out, parasites, hosts, hidden_cells=hidden_cells,
                target_flips=len(hidden_cells),
                host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
                flip_cost_matrix=flip_cost_matrix,
                tol=0, max_iter=20
            )
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
