# simulate.py
import argparse
import os
import numpy as np
from Bio import Phylo
import sys
import os

# Add parent directory to sys.path
sys.path.append(os.path.abspath(".."))

# Now import the module
from simulation_utils import *
import math
import json


def smooth_scores(scores, window=5, poly=2):
    if len(scores) < window:
        return scores
    return savgol_filter(scores, window_length=window, polyorder=poly)

def parsimony_vs_flips(out, parasites, hosts, host_W_matrices, par_W_matrices,
                       flip_cost_matrix, steps, lower, upper):
    """
    Compute total parsimony and number of flips across a range of λ values.
    Returns arrays of flips, scores, and lambdas.
    """
    lambdas = np.linspace(lower, upper, steps)
    scores = []
    flips_counts = []

    for lam in lambdas:
        cut_result = solve_network_cut(
            out, host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
            flip_cost_matrix=flip_cost_matrix,
            lambda_param=lam
        )

        # Count flips
        flips_counts.append(len(cut_result["flips"]))

        # total parsimony across parasite×host
        total = 0
        for p in parasites:
            leaf_states = {h: cut_result["new_cell_state"][(p,h)] for h in hosts}
            total += sankoff(out["host_trees"][0], leaf_states, np.array([[0,1],[1,0]]))
        for h in hosts:
            leaf_states = {p: cut_result["new_cell_state"][(p,h)] for p in parasites}
            total += sankoff(out["par_trees"][0], leaf_states, np.array([[0,1],[1,0]]))

        scores.append(total)

    return np.array(flips_counts), np.array(scores), np.array(lambdas)


def find_elbow_parsimony_flips(out, parasites, hosts, host_W_matrices,
                               par_W_matrices, flip_cost_matrix,
                               lower, upper, steps=50, outdir="experiments"):
    # Run sweep
    flips, scores, lambdas = parsimony_vs_flips(
        out, parasites, hosts, host_W_matrices, par_W_matrices,
        flip_cost_matrix, steps=steps, lower=lower, upper=upper
    )

    # Smooth scores
    smooth = smooth_scores(scores, window=3, poly=2)

    # Elbow detection in flip-space
    kneedle = KneeLocator(
        flips, smooth,
        curve="convex", direction="decreasing"
    )
    elbow_flip = kneedle.knee
    if elbow_flip is not None:
        # Find closest λ to the elbow flip
        idx = (np.abs(flips - elbow_flip)).argmin()
        elbow_lambda = lambdas[idx]
    else:
        # take lambda with max second derivative
        elbow_flip, elbow_lambda = None, lambdas[-1]

    print(f"[Elbow detection] flip ≈ {elbow_flip}, λ ≈ {elbow_lambda}")

    # Plot
    plt.figure(figsize=(10, 10))
    plt.plot(flips, scores, "bo-", label="Pareto Front")

    if elbow_flip is not None:
        # Find the corresponding y value for the elbow point
        idx = int(np.argmin(np.abs(np.array(flips) - elbow_flip)))
        elbow_y = scores[idx]
        # Mark the elbow with a red five-point star
        plt.plot(elbow_flip, elbow_y,
                 marker='*', markersize=22,
                 markerfacecolor='red', markeredgecolor='k',
                 label="CoEvoLink Solution")

    # Annotate points with λ
    # for f, s, lam in zip(flips, scores, lambdas):
    #     plt.text(f, s, f"{lam:.2f}", fontsize=11,
    #              ha="right", va="bottom", rotation=45)

    plt.xlabel("Prediction Cost", fontsize=24)
    plt.ylabel("Parsimony Score", fontsize=24)
    # incresae tick font size
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=24)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "parsimony_vs_flips_elbow.svg"))
    plt.close()

    return elbow_flip, elbow_lambda, flips, scores, lambdas



def read_nexus_reorder(nexus_file, parasite_tree_file, host_names):
    parasite_tree = Phylo.read(parasite_tree_file, "newick")
    parasite_order = [t.name for t in parasite_tree.get_terminals()]

    cell_state = {}
    hidden_cells = []  # store where the 1’s were

    # Parse NEXUS file manually
    reading = False
    rows = {}
    with open(nexus_file) as f:
        for line in f:
            line = line.strip()
            if line.lower().startswith("matrix"):
                reading = True
                continue
            if reading:
                if line.startswith(";"):
                    break
                if not line:
                    continue
                parts = line.split()
                taxon = parts[0]
                states = parts[1]
                rows[taxon] = states

    # Reorder according to parasite tree
    for parasite in parasite_order:
        if parasite not in rows:
            raise ValueError(f"Parasite {parasite} not in Nexus")
        states = rows[parasite]
        if len(states) != len(host_names):
            raise ValueError(f"Number of characters for {parasite} does not match number of hosts")

        for i, h in enumerate(host_names):
            s = int(states[i])
            if s == 1:  # treat 1’s as hidden positives
                cell_state[(parasite, h)] = 0  # mask to 0
                hidden_cells.append((parasite, h))
            elif s == 2:  # keep 2 as real positive
                cell_state[(parasite, h)] = 1
            else:  # 0 or others
                cell_state[(parasite, h)] = 0

    return cell_state, hidden_cells
def parse_braga_settings(settings_file):
    vals = {}
    with open(settings_file) as f:
        for line in f:
            if "," not in line or line.startswith("variable"):
                continue
            key, value = line.strip().split(",", 1)
            try:
                vals[key] = float(value)
            except ValueError:
                vals[key] = value.strip('"')
    return vals
def rates_from_settings(settings_file):
    vals = parse_braga_settings(settings_file)

    clock = float(vals["clock_host"])
    beta = float(vals["phy_scale[1]"])
    s01 = float(vals["switch_rate_0_to_1"])
    s12 = float(vals["switch_rate_1_to_2"])
    s10 = float(vals["switch_rate_1_to_0"])
    s21 = float(vals["switch_rate_2_to_1"])

    r01 = (s01 + s12) * math.exp(clock - beta)
    r10 = (s10 + s21) * math.exp(clock)

    return r01, r10

def matrix_builder(out, r01_h, r10_h, r01_p, r10_p, scale=1.0):
    # Build weight matrices for each parasite’s host tree
    host_W_matrices = []
    for t in out["host_trees"]:
        W, nodes = build_weight_matrix(t, r01_h, r10_h, scale=scale)
        host_W_matrices.append((W, nodes))

    # Build weight matrices for each host’s parasite tree
    par_W_matrices = []
    for t in out["par_trees"]:
        W, nodes = build_weight_matrix(t, r01_p, r10_p, scale=scale)
        par_W_matrices.append((W, nodes))

    # Flip cost matrix
    flip_cost_matrix = build_flip_cost_matrix(len(out["par_leaves"]), len(out["host_leaves"]), cost=1.0)
    return host_W_matrices, par_W_matrices, flip_cost_matrix





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--mode", choices=["ml", "default"], default="ml")
    parser.add_argument("--r01_p", type=float, default=0.5)
    parser.add_argument("--r10_p", type=float, default=0.5)
    parser.add_argument("--r01_h", type=float, default=0.5)
    parser.add_argument("--r10_h", type=float, default=0.5)
    parser.add_argument("--host_tree", type=str, required=True, help="Host tree file in Newick format")
    parser.add_argument("--virus_tree", type=str, required=True, help="Virus/parasite tree file in Newick format")
    parser.add_argument("--corrupt", type=int, default=0)
    parser.add_argument("--lambda_param", type=float, default=0.7, help="Lambda parameter for network cut")
    parser.add_argument("--outdir", type=str, default="experiments", help="Output directory for results")  # NEW
    parser.add_argument("--nexus_file", type=str, default=None,
                    help="Optional Nexus file containing parasite x host states")
    parser.add_argument("--original_metrics", type=str, default=None, help="File to save precision/recall/F1 metrics as JSON")
    parser.add_argument("--elbow_metrics", type=str, default=None, help="File to save precision/recall/F1 metrics as JSON")
    parser.add_argument("--braga", type=str, default=None, help="BRAGA settings file to compute rates")

    args = parser.parse_args()

    seed = args.seed
    r01_p = args.r01_p
    r10_p = args.r10_p
    r01_h = args.r01_h
    r10_h = args.r10_h
    mode = args.mode
    corrupt = args.corrupt
    outdir = args.outdir   # NEW

    os.makedirs(outdir, exist_ok=True)   # NEW: central output dir

    # If nexus_file is given, read it
    if args.nexus_file:
        # Read trees for structure only
        base_host_tree = Phylo.read(args.host_tree, "newick")
        base_par_tree  = Phylo.read(args.virus_tree, "newick")
        
        base_host_tree = rescale_tree(base_host_tree)
        base_par_tree  = rescale_tree(base_par_tree)
        
        host_leaves = [t.name for t in base_host_tree.get_terminals()]
        par_leaves  = [t.name for t in base_par_tree.get_terminals()]
        
        # Reorder & combine states from nexus
        cell_state, hidden_cells = read_nexus_reorder(args.nexus_file, args.virus_tree, host_leaves)
        
        out = {
            "host_trees": [copy.deepcopy(base_host_tree) for _ in range(len(par_leaves))],
            "par_trees":  [copy.deepcopy(base_par_tree)  for _ in range(len(host_leaves))],
            "cell_state": cell_state,
            "host_leaves": host_leaves,
            "par_leaves": par_leaves
        }

        host_leaves = out["host_leaves"]
        par_leaves = out["par_leaves"]

        mat, parasites, hosts = get_interaction_matrix(out)

    else:
        # Original simulation
        out = simulate_network(args.host_tree, args.virus_tree, r01_p, r10_p, r01_h, r10_h, seed, scale=100.0)

        host_leaves = out["host_leaves"]
        par_leaves = out["par_leaves"]

        score = total_parsimony(
            host_tree_file=args.host_tree,
            parasite_tree_file=args.virus_tree,
            cell_state=out["cell_state"]
        )
        mat, parasites, hosts = get_interaction_matrix(out)

        # if mode == "ml":
        #     max_ll = 0.0
        #     for i, host_tree in enumerate(out["host_trees"]):
        #         for p in par_leaves:
        #             leaf_states = {h: out["cell_state"][(p,h)] for h in host_leaves}
        #             max_ll += sankoff_ml_loglik(host_tree, leaf_states, r01_h, r10_h)
        #     for j, par_tree in enumerate(out["par_trees"]):
        #         for h in host_leaves:
        #             leaf_states = {p: out["cell_state"][(p,h)] for p in par_leaves}
        #             max_ll += sankoff_ml_loglik(par_tree, leaf_states, r01_p, r10_p )

        #     current_ll = total_current_ll(
        #         out["host_trees"], out["par_trees"], out["cell_state"], host_leaves, par_leaves
        #     )
        #     ratio = current_ll / max_ll

        #     if 1:
        #         os.makedirs("ml_parsimony_matrices", exist_ok=True)
        #         filename = os.path.join(
        #             "ml_parsimony_matrices", 
        #             f"seed{seed}_ratio{ratio:.4f}.png"
        #         )
        #         plot_matrix(mat, parasites, hosts, filename=filename)

        # else:
        #     os.makedirs("default_parsimony_matrices", exist_ok=True)
        #     filename = os.path.join(
        #         "default_parsimony_matrices", 
        #         f"seed{seed}_{score}.png"
        #     )
        #     plot_matrix(mat, parasites, hosts, filename=filename)
    

    # Save original and flipped matrices as images
    # os.makedirs("experiments", exist_ok=True)

    # Original
    # plot_matrix(mat, parasites, hosts,filename=os.path.join(outdir, "original_matrix.png"))

    if corrupt == 0:
        lambda_param = args.lambda_param
    
    flip_cost = 1

    host_W_matrices, par_W_matrices, flip_cost_matrix = matrix_builder(out, r01_h, r10_h, r01_p, r10_p, scale=100.0)
        
    # ---------------------------
    # Step 1: Create corrupted matrix (random 1→0 flips)
    # ---------------------------


    if corrupt > 0:
        random_flips = min(corrupt, np.sum(mat) - 1)
        corrupt_mat, hidden = corrupt_matrix(mat, num_flips=random_flips, seed=seed)

        # Save hidden cells for metrics
        hidden_cells = [(parasites[i], hosts[j]) for i, j in hidden]

        highlight_corrupted = {"corrupted": hidden_cells}
        # plot_matrix(corrupt_mat, parasites, hosts,
        #             filename=os.path.join(outdir, "corrupted.png"),
        #             highlight=highlight_corrupted)

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
        if args.braga:
            if args.nexus_file is None or args.braga is False:
                raise ValueError("BRAGA mode requires --nexus_file and --braga settings file")

            # Compute rates
            r01, r10 = rates_from_settings(args.braga)
            r01_p = r01_h = r01
            r10_p = r10_h = r10
            # print (r01, r10)

            if np.sum(mat) == 0:
        # force one positive
                i, j = random.randrange(mat.shape[0]), random.randrange(mat.shape[1])
                mat[i, j] = 1
                parasite, host = parasites[i], hosts[j]
                out["cell_state"][(parasite, host)] = 1
                print(f"[WARN] Matrix was all zeros. Forced a 1 at ({parasite}, {host}).")


            highlight_hidden = {"corrupted": hidden_cells}
            # plot_matrix(mat, parasites, hosts,filename=os.path.join(outdir, "corrupted.png"),
            # highlight=highlight_hidden)
            lambda_param, cut_result = binary_search_lambda(out, parasites, hosts, hidden_cells=hidden_cells,
            target_flips=len(hidden_cells), host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
            flip_cost_matrix=flip_cost_matrix, tol=0, max_iter=20)
            # cut_result = solve_network_cut(out, lambda_param=lambda_param, flip_cost=flip_cost, r01_p=r01_p, r10_p=r10_p, r01_h=r01_h, r10_h=r10_h)

        else:
            highlight_hidden = {"corrupted": hidden_cells}
            # plot_matrix(mat, parasites, hosts,filename=os.path.join(outdir, "corrupted.png"),
            # highlight=highlight_hidden)
            cut_result = solve_network_cut(out, host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
                                            flip_cost_matrix=flip_cost_matrix, lambda_param=lambda_param)
    # ---------------------------
    # Step 2: Run network cut recovery on corrupted input
    # ---------------------------
    

    new_mat = cut_result["new_matrix"]
    new_cell_state = cut_result["new_cell_state"]
    flips = cut_result["flips"]

    highlight_flipped = {"flipped": [(p, h) for p, h, old, new in flips]}
    print(f"Number of flips performed: {len(flips)}")
    # for p, h, old, new in flips:
    #     print(f"Cell ({p},{h}): {old} -> {new}")

    # ---------------------------
    # Step 3: Save recovered matrix with algorithm flips
    # ---------------------------
    # plot_matrix(new_mat, parasites, hosts,
    #             filename=os.path.join(outdir, "flipped_matrix.svg"),
    #             highlight=highlight_flipped)






    # inside main(), after you have flips and hidden_cells
    metrics = compute_metrics(hidden_cells, flips)
    print("Metrics:", metrics)


    elbow_flip, elbow_lambda, flips, scores, lambdas = find_elbow_parsimony_flips(
        out, parasites, hosts, host_W_matrices, par_W_matrices, flip_cost_matrix, lower=0.0, upper=1.0, steps=50, outdir=outdir
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
    # plot_matrix(cut_result_elbow["new_matrix"], parasites, hosts,
    #             filename=os.path.join(outdir, "flipped_matrix_elbow.png"),
    #             highlight=highlight_flipped)

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


# for paper simulation figure I did seed 31 r 0.1, 0.8 corrupt 200
