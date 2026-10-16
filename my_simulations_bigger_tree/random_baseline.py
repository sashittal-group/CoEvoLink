# random_baseline.py
#
# Simple baselines for the recovery experiments (Reviewer 3, comment 7).
#
# The simulated data are generated exactly as in simulation.py -- same trees,
# same seeds, same transition rates, same corruption -- so the true interaction
# matrix and the set of hidden (1 -> 0) cells are bit-for-bit identical to the
# ones CoEvoLink was evaluated on. Only the recovery step changes: instead of
# solving the network cut, we pick cells among the 0s of the corrupted matrix
# and flip them to 1, using
#
#   random      : cells chosen uniformly at random;
#   greedy      : cells belonging to the hosts with the largest observed number
#                 of interactions first, i.e. missing entries are attributed to
#                 the most connected hosts;
#   greedy_prod : cells ranked by the product of the observed degrees of their
#                 virus and their host, so a missing entry is attributed to a
#                 cell only when both of its margins are well connected.
#
# greedy_prod is the degree heuristic suggested in review.  It scores every 0
# cell once on the observed matrix and takes the top-k, rather than flipping
# the single best cell and re-scoring: re-scoring is self-reinforcing, because
# a flip raises the score of every other cell in the same row and column, and
# it measures worse.
#
# All baselines use only the observed (corrupted) matrix -- no phylogeny, no
# genomic scores.
#
# Two budgets, matching the two CoEvoLink settings reported in the paper:
#   supervised   : the number of missing entries is known, so the baseline gets
#                  the same budget CoEvoLink gets through binary_search_lambda
#                  (target_flips = |hidden cells|).
#   unsupervised : the budget is not known, so the baseline gets the number of
#                  flips the elbow (Pareto) solution actually made, read from
#                  the corresponding CoEvoLink elbow_metrics.json (TP + FP).
#
# One invocation handles every corruption level for a given (seed, rate) pair,
# because simulate_network() does not depend on the corruption level and is by
# far the most expensive step.

import argparse
import json
import os
import sys

import numpy as np

sys.path.append(os.path.abspath(".."))

from simulation_utils import corrupt_matrix, get_interaction_matrix, simulate_network


# Method order fixes the RNG tag of the two original methods, so it is not
# cosmetic: "random" and "greedy" must stay in slots 0 and 1.
METHODS = ["random", "greedy", "greedy_prod"]

# base tag for methods added after the original run, kept clear of tags 0-3
RNG_TAG_NEW = 1000


def metrics_from_counts(TP, FP, FN):
    """Same definitions as simulation_utils.compute_metrics, on counts."""
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return {"precision": precision, "recall": recall, "f1": f1,
            "TP": TP, "FP": FP, "FN": FN}


def score_and_average(pick_fn, hidden_mask_flat, budget, n_zeros, repeats):
    """Average precision/recall/F1 over `repeats` independent picks of `budget` cells."""
    n_hidden = int(hidden_mask_flat.sum())
    budget = int(min(budget, n_zeros))
    if budget == 0:
        return metrics_from_counts(0, 0, n_hidden), budget

    acc = {"precision": 0.0, "recall": 0.0, "f1": 0.0, "TP": 0.0, "FP": 0.0, "FN": 0.0}
    for _ in range(repeats):
        picked = pick_fn(budget)
        TP = int(hidden_mask_flat[picked].sum())
        m = metrics_from_counts(TP, budget - TP, n_hidden - TP)
        for k in acc:
            acc[k] += m[k]
    return {k: v / repeats for k, v in acc.items()}, budget


def make_random_picker(zero_idx, rng):
    """Uniformly random cells among the 0s of the corrupted matrix."""
    def pick(budget):
        return rng.choice(zero_idx, size=budget, replace=False)
    return pick


def make_greedy_picker(zero_idx, zero_scores, rng):
    """
    Top-scoring 0 cells, ties broken at random.

    The caller supplies the score of each candidate cell; see the two score
    vectors built in main().  Degrees are integers, so scores tie constantly
    and the random permutation below is what keeps the tie-breaking unbiased.
    """
    def pick(budget):
        order = rng.permutation(len(zero_idx))          # random tie-breaking
        ranked = order[np.argsort(-zero_scores[order], kind="stable")]
        return zero_idx[ranked[:budget]]
    return pick


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--host_tree", type=str, required=True)
    parser.add_argument("--virus_tree", type=str, required=True)
    parser.add_argument("--corrupts", type=int, nargs="+", required=True,
                        help="Corruption levels to evaluate for this (seed, rate) pair")
    parser.add_argument("--r01_p", type=float, default=0.5)
    parser.add_argument("--r10_p", type=float, default=0.5)
    parser.add_argument("--r01_h", type=float, default=0.5)
    parser.add_argument("--r10_h", type=float, default=0.5)
    parser.add_argument("--outroot", type=str, default="results_baseline",
                        help="Root of the baseline result tree")
    parser.add_argument("--coevolink_root", type=str, default="results/consideringGenomic",
                        help="Root of the CoEvoLink results, used to read the elbow flip budget")
    parser.add_argument("--repeats", type=int, default=5,
                        help="Independent draws averaged per configuration")
    parser.add_argument("--baseline_seed_offset", type=int, default=100000,
                        help="Offset keeping the baseline RNG independent of the simulation RNG")

    args = parser.parse_args()

    # Path fragment shared by both result trees, e.g. r010.1_r100.1
    rate_dir = f"r01{args.r01_p}_r10{args.r10_p}"

    # ------------------------------------------------------------------
    # Step 1: regenerate the exact same ground-truth matrix as simulation.py
    # ------------------------------------------------------------------
    out = simulate_network(args.host_tree, args.virus_tree,
                           args.r01_p, args.r10_p, args.r01_h, args.r10_h,
                           args.seed, scale=100.0)
    mat, parasites, hosts = get_interaction_matrix(out)
    n_ones = int(mat.sum())
    P, H = mat.shape

    for corrupt in args.corrupts:
        outdir = os.path.join(args.outroot, str(args.seed), f"corrupt{corrupt}", rate_dir)
        os.makedirs(outdir, exist_ok=True)

        # --------------------------------------------------------------
        # Step 2: identical corruption (same seed -> same hidden cells)
        # --------------------------------------------------------------
        random_flips = min(corrupt, n_ones - 1)
        corrupt_mat, hidden = corrupt_matrix(mat, num_flips=random_flips, seed=args.seed)

        hidden_mask = np.zeros((P, H), dtype=bool)
        for i, j in hidden:
            hidden_mask[i, j] = True
        hidden_mask_flat = hidden_mask.ravel()
        n_hidden = len(hidden)

        zero_idx = np.flatnonzero(corrupt_mat.ravel() == 0)
        n_zeros = len(zero_idx)

        # Greedy scores, both computed on the corrupted matrix only:
        #   greedy      -> observed degree of the host
        #   greedy_prod -> observed virus degree x observed host degree
        host_degree = corrupt_mat.sum(axis=0)
        virus_degree = corrupt_mat.sum(axis=1)
        greedy_scores = {
            "greedy": np.tile(host_degree, P)[zero_idx].astype(float),
            "greedy_prod": np.outer(virus_degree,
                                    host_degree).ravel()[zero_idx].astype(float),
        }

        # --------------------------------------------------------------
        # Step 3: budgets -- known number of missing entries (supervised),
        #         number of flips made by the CoEvoLink elbow (unsupervised)
        # --------------------------------------------------------------
        elbow_path = os.path.join(args.coevolink_root, str(args.seed),
                                  f"corrupt{corrupt}", rate_dir, "elbow_metrics.json")
        if os.path.exists(elbow_path):
            with open(elbow_path) as f:
                elbow = json.load(f)
            elbow_budget = int(elbow["TP"] + elbow["FP"])
            elbow_source = "coevolink_elbow_flips"
        else:
            # No CoEvoLink elbow run to match: fall back to the supervised budget
            # and record that fact so it can be filtered out downstream.
            elbow_budget = n_hidden
            elbow_source = "missing_elbow_fallback_n_hidden"
            print(f"[WARN] no elbow metrics at {elbow_path}; falling back to n_hidden")

        budgets = {
            "supervised":   (n_hidden,     "n_hidden"),
            "unsupervised": (elbow_budget, elbow_source),
        }

        common = {
            "seed": args.seed,
            "corrupt": corrupt,
            "r01": args.r01_p,
            "r10": args.r10_p,
            "n_ones": n_ones,
            "n_hidden": n_hidden,
            "n_zeros": n_zeros,
            "missing_rate": n_hidden / n_ones if n_ones else 0.0,
            "repeats": args.repeats,
        }

        # A distinct, reproducible stream per (config, mode, method)
        def rng_for(tag):
            return np.random.default_rng(
                [args.baseline_seed_offset, args.seed, corrupt,
                 int(1000 * args.r01_p), int(1000 * args.r10_p), tag]
            )

        line = [f"seed={args.seed} corrupt={corrupt} r={args.r01_p}/{args.r10_p} "
                f"ones={n_ones} hidden={n_hidden} zeros={n_zeros}"]

        for mode_tag, (mode, (budget, budget_source)) in enumerate(budgets.items()):
            for method_tag, method in enumerate(METHODS):
                # tags 0-3 belong to the two original methods; greedy_prod draws
                # from a separate stream so their numbers are reproduced exactly
                tag = (RNG_TAG_NEW + mode_tag if method == "greedy_prod"
                       else 2 * mode_tag + method_tag)
                rng = rng_for(tag)
                pick_fn = (make_random_picker(zero_idx, rng) if method == "random"
                           else make_greedy_picker(zero_idx,
                                                   greedy_scores[method], rng))
                res, used = score_and_average(pick_fn, hidden_mask_flat,
                                              budget, n_zeros, args.repeats)
                res.update(common)
                res["budget"] = used
                res["budget_source"] = budget_source
                res["method"] = method
                with open(os.path.join(outdir, f"{mode}_{method}_metrics.json"), "w") as f:
                    json.dump(res, f, indent=2)
                line.append(f"{mode}/{method} F1={res['f1']:.4f} (k={used})")

        print(" | ".join(line))


if __name__ == "__main__":
    main()
