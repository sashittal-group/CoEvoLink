#!/usr/bin/env python
"""
parametric_pareto.py -- enumerate the CoEvoLink Pareto front by parametric
maximum flow (Module 1; Theorem 2 and Algorithm 2 of the Supplemental Material).

Balanced PIP at weight lambda is the minimum s-t cut of G(A', T_v, T_h), whose
tree edges have capacity lambda * w and whose source edges have capacity
(1 - lambda) * c (see simulation_utils.solve_network_cut). Dividing every
capacity by lambda and writing mu = (1 - lambda) / lambda, each capacity becomes
affine in mu,

    cap(mu) = a + b * mu,

with b = 0 on tree edges and b = c_ij on source edges. The value of a fixed cut
is therefore the line  S_w(A) + mu * c(A, A'),  and the optimum over all cuts is
the lower envelope of these lines. The pieces of that envelope are exactly the
Pareto optimal interaction matrices on the lower convex hull of the front.

They are found with the divide-and-conquer of Gallo, Grigoriadis and Tarjan
(1989): between two known solutions the only candidate breakpoint is the value
mu* at which their two lines intersect, so one minimum cut at mu* either finds
a new matrix (recurse on both sides) or certifies mu* as a breakpoint. Because
only the source capacities depend on mu, and they grow with it, the network is
monotone and its minimum cuts are nested, so each recursive call only has to
solve the vertices lying strictly between the two known cuts, with the rest
contracted into s and t.

Implementation notes
--------------------
* Minimum cuts are computed with networkx (preflow-push), the same routine
  solve_network_cut uses, and with the same tie-breaking: networkx returns the
  largest source side among the minimum cuts. The matrix this script reports
  for an interval of lambda is therefore the matrix solve_network_cut returns
  for any lambda inside that interval (check it with --check_grid).
* Each contracted subnetwork is solved from scratch rather than warm-started
  from the parent flow, so this reference implementation does not attain the
  single-maximum-flow running time of Theorem 2; it returns the same matrices
  and breakpoints.
* Elbow (Module 2): kneedle depends on which points of the front it is given.
  The results in the paper apply it to the solutions at 50 evenly spaced
  values of lambda (find_elbow_parsimony_flips), and --elbow_grid 50, the
  default, reproduces that choice from the enumerated front without solving
  any further minimum cuts. --elbow_grid 0 applies kneedle to every matrix on
  the front instead.
* Cells whose flip cost is exactly 100.0 are handled as in solve_network_cut,
  with an edge to the sink of constant capacity 100 (the VIRION scripts use it
  to pin known hosts). After dividing by lambda that capacity is 100 * (1 + mu):
  the value of a cut is still affine in mu, but the network is no longer
  monotone. The script detects this and computes every minimum cut on the full
  network instead of a contracted one; the envelope recursion is still exact.

Examples
--------
Regenerate one instance of the Mk-model simulation (the same inputs as
my_simulations_bigger_tree/simulation.py) and enumerate its Pareto front:

    python parametric_pareto.py --simulate --seed 42 --corrupt 1000 \\
        --r01 0.1 --r10 0.4 \\
        --host_tree my_simulations_bigger_tree/hosts_tree_big.newick \\
        --virus_tree my_simulations_bigger_tree/viruses_tree_big.newick \\
        --outdir parametric_out/seed42

Run on your own data (rows = viruses, columns = hosts, entries 0/1):

    python parametric_pareto.py --matrix A.csv --costs C.csv \\
        --host_tree hosts.newick --virus_tree viruses.newick --outdir out
"""

import argparse
import bisect
import copy
import json
import os
import sys

import numpy as np
import networkx as nx
import pandas as pd
from networkx.algorithms.flow import preflow_push
from Bio import Phylo

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simulation_utils import (  # noqa: E402
    build_flip_cost_matrix,
    build_flip_cost_matrix_considering_genomic,
    build_weight_matrix,
    corrupt_matrix,
    get_interaction_matrix,
    rescale_tree,
    sankoff,
    simulate_network,
    solve_network_cut,
)

SOURCE, SINK = 0, 1
PIN_COST = 100.0  # sentinel flip cost treated specially by solve_network_cut
INF = float("inf")


# ---------------------------------------------------------------------------
# Parametric network G(mu)
# ---------------------------------------------------------------------------
class ParametricNetwork:
    """G(A', T_v, T_h) with capacities a + b * mu, mu = (1 - lambda) / lambda.

    Takes exactly the inputs of simulation_utils.solve_network_cut and builds
    the same graph, edge for edge.
    """

    def __init__(self, out, host_W_matrices, par_W_matrices, flip_cost_matrix):
        self.out = out
        self.hosts = list(out["host_leaves"])
        self.parasites = list(out["par_leaves"])
        self.cost = np.asarray(flip_cost_matrix, dtype=float)
        self.cell_state = out["cell_state"]

        self._ids = {}
        self._edges = {}  # (u, v) -> (a, b); later additions overwrite, as in networkx
        self._node_id("SOURCE")
        self._node_id("SINK")

        for p_idx, (p, t) in enumerate(zip(self.parasites, out["host_trees"])):
            W, nodes = host_W_matrices[p_idx]
            self._add_tree(t, W, nodes, lambda leaf: f"CELL_{p}_{leaf}")
        for h_idx, (h, t) in enumerate(zip(self.hosts, out["par_trees"])):
            W, nodes = par_W_matrices[h_idx]
            self._add_tree(t, W, nodes, lambda leaf: f"CELL_{leaf}_{h}")

        self.pinned = np.zeros(self.cost.shape, dtype=bool)
        self.observed = np.zeros(self.cost.shape, dtype=int)
        self.cell_ids = np.zeros(self.cost.shape, dtype=np.int64)
        for i, p in enumerate(self.parasites):
            for j, h in enumerate(self.hosts):
                node = f"CELL_{p}_{h}"
                self.cell_ids[i, j] = self._node_id(node)
                state = self.cell_state[(p, h)]
                self.observed[i, j] = state
                c = self.cost[i, j]
                if state == 0:
                    if c == PIN_COST:
                        # constant capacity c in lambda-space -> c * (1 + mu)
                        self.pinned[i, j] = True
                        self._add(node, "SINK", c, c)
                        self._add("SINK", node, c, c)
                    else:
                        self._add("SOURCE", node, 0.0, c)
                else:
                    self._add(node, "SINK", INF, 0.0)
                    self._add("SINK", node, INF, 0.0)

        self.par_root_ids = [self._node_id(t.root) for t in out["par_trees"]]

        self.n_nodes = len(self._ids)
        keys = list(self._edges.keys())
        self.eu = np.array([u for u, _ in keys], dtype=np.int64)
        self.ev = np.array([v for _, v in keys], dtype=np.int64)
        ab = np.array(list(self._edges.values()), dtype=float)
        self.ea, self.eb = ab[:, 0], ab[:, 1]
        self.tree_edge = (self.eu > SINK) & (self.ev > SINK)

        # Monotone iff mu only enters source capacities, and only increasingly.
        param = self.eb != 0
        ok = ((self.eu == SOURCE) & (self.eb > 0)) | ((self.ev == SINK) & (self.eb < 0))
        self.monotone = bool(np.all(ok[param]))

        self.total_tree_weight = float(self.ea[self.tree_edge].sum()) / 2.0
        positive = self.cost[(self.observed == 0) & (self.cost > 0)]
        self.min_cost = float(positive.min()) if positive.size else 1.0

    def _node_id(self, key):
        if key not in self._ids:
            self._ids[key] = len(self._ids)
        return self._ids[key]

    def _add(self, u, v, a, b):
        self._edges[(self._node_id(u), self._node_id(v))] = (a, b)

    def _add_tree(self, tree, W, nodes, cell_name):
        index = {id(n): k for k, n in enumerate(nodes)}
        for parent in tree.find_clades(order="preorder"):
            for child in parent.clades:
                w = W[index[id(parent)], index[id(child)]]
                node = cell_name(child.name) if child.is_terminal() else child
                self._add(parent, node, w, 0.0)
                self._add(node, parent, w, 0.0)

    # -- cuts ---------------------------------------------------------------
    def line(self, in_R):
        """(intercept, slope) of the cut with source side in_R, as a function of mu."""
        crossing = in_R[self.eu] & ~in_R[self.ev]
        return float(self.ea[crossing].sum()), float(self.eb[crossing].sum())

    def min_cut(self, mu, in_S=None, in_T=None):
        """Largest-source-side minimum cut of G(mu) with in_S fixed on the
        source side and in_T fixed on the sink side (contracted into s and t)."""
        n = self.n_nodes
        if in_S is None:
            in_S = np.zeros(n, dtype=bool)
            in_S[SOURCE] = True
        if in_T is None:
            in_T = np.zeros(n, dtype=bool)
            in_T[SINK] = True
        mid = ~in_S & ~in_T
        if not mid.any():
            return in_S.copy()

        cap = self.ea + mu * self.eb
        mu_mid, mv_mid = mid[self.eu], mid[self.ev]
        s, t = -1, -2

        G = nx.DiGraph()
        G.add_node(s)
        G.add_node(t)
        G.add_nodes_from(np.flatnonzero(mid).tolist())

        sel = mu_mid & mv_mid & (cap > 0)
        G.add_edges_from(
            (u, v, {"capacity": c})
            for u, v, c in zip(self.eu[sel].tolist(), self.ev[sel].tolist(), cap[sel].tolist())
        )
        for side, sel, endpoint in (
            (s, in_S[self.eu] & mv_mid, self.ev),
            (t, mu_mid & in_T[self.ev], self.eu),
        ):
            agg = np.zeros(n)
            np.add.at(agg, endpoint[sel], cap[sel])
            for x in np.flatnonzero(agg > 0).tolist():
                if side == s:
                    G.add_edge(s, x, capacity=agg[x])
                else:
                    G.add_edge(x, t, capacity=agg[x])

        # Same partition rule as nx.minimum_cut (the sink side is every vertex
        # that can still reach t in the residual network), except that an arc
        # counts as saturated when its residual capacity is within a relative
        # tolerance of zero: nx.minimum_cut tests flow == capacity exactly,
        # which on float capacities can leave a saturated arc in the residual
        # network and return a cut that is not minimum.
        R = preflow_push(G, s, t, value_only=True)
        finite = [c for _, _, c in G.edges(data="capacity") if c < INF]
        eps = 1e-9 * max([1.0] + finite)
        sink_side, frontier = {t}, [t]
        while frontier:
            v = frontier.pop()
            for u in R.predecessors(v):
                if u not in sink_side and R[u][v]["capacity"] - R[u][v]["flow"] > eps:
                    sink_side.add(u)
                    frontier.append(u)
        in_R = in_S.copy()
        in_R[[x for x in G if x >= 0 and x not in sink_side]] = True
        return in_R

    # -- reading a cut ------------------------------------------------------
    def solution(self, in_R):
        """Same fields as solve_network_cut, plus the two objectives."""
        new_matrix = np.where(in_R[self.cell_ids], 0, 1).astype(int)
        new_cell_state, flips = {}, []
        for i, p in enumerate(self.parasites):
            for j, h in enumerate(self.hosts):
                new, old = int(new_matrix[i, j]), int(self.observed[i, j])
                new_cell_state[(p, h)] = new
                if new != old:
                    flips.append((p, h, old, new))
        root_states = {
            h: (0 if in_R[r] else 1) for h, r in zip(self.hosts, self.par_root_ids)
        }
        crossing = in_R[self.eu] & ~in_R[self.ev]
        flipped = (new_matrix == 1) & (self.observed == 0)
        return {
            "flips": flips,
            "new_matrix": new_matrix,
            "new_cell_state": new_cell_state,
            "root_states": root_states,
            "weighted_parsimony": float(self.ea[crossing & self.tree_edge].sum()),
            "prediction_cost": float(self.cost[flipped & ~self.pinned].sum()),
        }


# ---------------------------------------------------------------------------
# Envelope enumeration (Algorithm 2 of the Supplemental Material)
# ---------------------------------------------------------------------------
def _intersect(l1, l2):
    return (l2[0] - l1[0]) / (l1[1] - l2[1])


def _lower_envelope(lines, tol):
    """Indices of the lines forming the lower envelope on mu >= 0, in order of
    increasing mu (decreasing slope), and the breakpoints between them."""
    order = sorted(range(len(lines)), key=lambda k: (-lines[k][1], lines[k][0]))
    dedup = []
    for k in order:  # equal slopes: keep the lowest intercept
        if dedup and abs(lines[dedup[-1]][1] - lines[k][1]) <= tol:
            continue
        dedup.append(k)
    hull = []
    for k in dedup:
        while len(hull) >= 2:
            l1, l2, l3 = lines[hull[-2]], lines[hull[-1]], lines[k]
            if _intersect(l1, l3) <= _intersect(l1, l2) + tol:
                hull.pop()
            else:
                break
        hull.append(k)
    while len(hull) >= 2 and _intersect(lines[hull[0]], lines[hull[1]]) <= tol:
        hull.pop(0)
    breaks = [_intersect(lines[a], lines[b]) for a, b in zip(hull, hull[1:])]
    return hull, breaks


def enumerate_pareto_parametric(out, host_W_matrices, par_W_matrices,
                                flip_cost_matrix, rel_tol=1e-9, verbose=False):
    """Enumerate the Pareto optimal interaction matrices on the lower convex
    hull of the front, with the values of lambda at which each one is optimal.

    Inputs are those of simulation_utils.solve_network_cut. Returns a list of
    dicts, ordered by increasing lambda (increasing number of flips); each has
    the fields of solve_network_cut plus lambda_low/lambda_high,
    mu_low/mu_high, weighted_parsimony and prediction_cost.
    """
    net = ParametricNetwork(out, host_W_matrices, par_W_matrices, flip_cost_matrix)
    contract = net.monotone
    if verbose:
        print(f"[parametric] {net.n_nodes} vertices, {len(net.eu)} arcs, "
              f"monotone={net.monotone}")

    mu_lo = 0.0
    mu_hi = 1.0 + net.total_tree_weight / min(net.min_cost, PIN_COST)

    cuts, lines, seen = [], [], set()

    def record(in_R):
        seen.add(in_R.tobytes())
        cuts.append(in_R)
        lines.append(net.line(in_R))
        return len(cuts) - 1

    def solve(mu, i1, i2):
        if contract:
            return net.min_cut(mu, in_S=cuts[i1], in_T=~cuts[i2])
        return net.min_cut(mu)

    k_lo = record(net.min_cut(mu_lo))
    k_hi = record(net.min_cut(mu_hi))

    stack = [(k_lo, k_hi)]
    while stack:
        i1, i2 = stack.pop()
        (a1, b1), (a2, b2) = lines[i1], lines[i2]
        if np.array_equal(cuts[i1], cuts[i2]) or abs(b1 - b2) <= rel_tol * max(1.0, abs(b1)):
            continue
        mu_star = (a2 - a1) / (b1 - b2)  # the only candidate breakpoint
        in_R = solve(mu_star, i1, i2)
        a, b = net.line(in_R)
        value, bound = a + mu_star * b, a1 + mu_star * b1
        if verbose:
            print(f"[parametric] mu*={mu_star:.6g}  cut={value:.6g}  line={bound:.6g}")
        if value >= bound - rel_tol * max(1.0, abs(bound)) or in_R.tobytes() in seen:
            continue  # mu_star is a breakpoint between consecutive solutions
        k = record(in_R)
        stack.append((i1, k))
        stack.append((k, i2))

    scale = max(1.0, max(abs(a) + abs(b) for a, b in lines))
    hull, breaks = _lower_envelope(lines, rel_tol * scale)

    # Interval of mu on which each solution is optimal (mu increasing).
    bounds = [mu_lo] + breaks + [INF]
    solutions = []
    for pos, k in enumerate(hull):
        m_lo, m_hi = bounds[pos], bounds[pos + 1]
        in_R = cuts[k]
        if not contract:
            # Without nesting, re-solve inside the interval so that ties are
            # broken exactly as solve_network_cut breaks them there.
            m_mid = (m_lo + m_hi) / 2 if m_hi < INF else max(mu_hi, 2 * m_lo + 1)
            in_R = net.min_cut(m_mid)
        sol = net.solution(in_R)
        sol["cut_intercept"], sol["cut_slope"] = net.line(in_R)
        sol.update({
            "mu_low": m_lo,
            "mu_high": m_hi,
            "lambda_low": 0.0 if m_hi == INF else 1.0 / (1.0 + m_hi),
            "lambda_high": 1.0 / (1.0 + m_lo),
        })
        solutions.append(sol)

    solutions.reverse()  # increasing lambda
    return solutions


def objective(solution, lam):
    """lambda * S_w + (1 - lambda) * c of a solution: the cut value that
    solve_network_cut reports at this lambda."""
    return lam * solution["cut_intercept"] + (1.0 - lam) * solution["cut_slope"]


def solution_at_lambda(solutions, lam):
    """The enumerated solution that solve_network_cut returns at this lambda."""
    if lam <= 0:
        return solutions[0]
    mu = (1.0 - lam) / lam
    # solutions are ordered by decreasing mu; at a breakpoint networkx keeps
    # the larger source side, i.e. the solution on the larger-mu side.
    mus = [s["mu_low"] for s in reversed(solutions)]  # increasing
    k = bisect.bisect_right(mus, mu) - 1
    return solutions[len(solutions) - 1 - k]


# ---------------------------------------------------------------------------
# Elbow (Module 2), as in find_elbow_parsimony_flips
# ---------------------------------------------------------------------------
def unweighted_parsimony(out, cell_state):
    """Parsimony score used for the Pareto-front plots in the simulation scripts."""
    C = np.array([[0, 1], [1, 0]])
    parasites, hosts = out["par_leaves"], out["host_leaves"]
    total = 0
    for p in parasites:
        total += sankoff(out["host_trees"][0], {h: cell_state[(p, h)] for h in hosts}, C)
    for h in hosts:
        total += sankoff(out["par_trees"][0], {p: cell_state[(p, h)] for p in parasites}, C)
    return total


def find_elbow(points):
    """Kneedle elbow on (flips, parsimony) points, as in the simulation scripts.
    Returns the index of the elbow point, or the last point if none is found."""
    from kneed import KneeLocator
    from scipy.signal import savgol_filter

    flips = np.array([f for f, _ in points])
    scores = np.array([s for _, s in points], dtype=float)
    smooth = savgol_filter(scores, window_length=3, polyorder=2) if len(scores) >= 3 else scores
    knee = KneeLocator(flips, smooth, curve="convex", direction="decreasing").knee
    if knee is None:
        return len(points) - 1
    return int(np.argmin(np.abs(flips - knee)))


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def build_out(host_tree_file, virus_tree_file, matrix):
    """The `out` structure of simulation_utils from two trees and a 0/1 matrix
    (pandas DataFrame, rows = viruses, columns = hosts)."""
    host_tree = rescale_tree(Phylo.read(host_tree_file, "newick"))
    par_tree = rescale_tree(Phylo.read(virus_tree_file, "newick"))
    hosts = [c.name for c in host_tree.get_terminals()]
    parasites = [c.name for c in par_tree.get_terminals()]
    matrix = matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    missing = (set(parasites) - set(matrix.index)) | (set(hosts) - set(matrix.columns))
    if missing:
        raise ValueError(f"tree leaves missing from the matrix: {sorted(missing)[:5]} ...")
    cell_state = {(p, h): int(matrix.at[p, h]) for p in parasites for h in hosts}
    return {
        "host_trees": [copy.deepcopy(host_tree) for _ in parasites],
        "par_trees": [copy.deepcopy(par_tree) for _ in hosts],
        "cell_state": cell_state,
        "host_leaves": hosts,
        "par_leaves": parasites,
    }


def weight_matrices(out, r01_h, r10_h, r01_p, r10_p, scale):
    host_W = [build_weight_matrix(t, r01_h, r10_h, scale=scale) for t in out["host_trees"]]
    par_W = [build_weight_matrix(t, r01_p, r10_p, scale=scale) for t in out["par_trees"]]
    return host_W, par_W


def simulated_instance(args):
    """One corrupted instance of the Mk-model simulation, built exactly as in
    my_simulations_bigger_tree/simulation.py."""
    out = simulate_network(args.host_tree, args.virus_tree, args.r01, args.r10,
                           args.r01, args.r10, args.seed, scale=args.scale)
    mat, parasites, hosts = get_interaction_matrix(out)
    host_W, par_W = weight_matrices(out, args.r01, args.r10, args.r01, args.r10, args.scale)
    costs = build_flip_cost_matrix(len(parasites), len(hosts), cost=1.0)
    hidden_cells = []
    if args.corrupt > 0:
        corrupt_mat, hidden = corrupt_matrix(mat, num_flips=min(args.corrupt, np.sum(mat) - 1),
                                             seed=args.seed)
        hidden_cells = [(parasites[i], hosts[j]) for i, j in hidden]
        if args.genomic_costs:
            costs = build_flip_cost_matrix_considering_genomic(
                parasites, hosts, hidden_cells, base_cost=1.0, hidden_mean=0.2,
                hidden_std=0.05, clip_min=0.01, clip_max=1.0, seed=args.seed)
        out["cell_state"] = {(p, h): int(corrupt_mat[i, j])
                             for i, p in enumerate(parasites) for j, h in enumerate(hosts)}
    return out, host_W, par_W, costs, hidden_cells


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Enumerate the CoEvoLink Pareto front by parametric maximum flow.")
    ap.add_argument("--host_tree", required=True, help="host tree (Newick)")
    ap.add_argument("--virus_tree", required=True, help="virus tree (Newick)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--matrix", help="observed 0/1 matrix CSV, rows = viruses, columns = hosts")
    src.add_argument("--simulate", action="store_true",
                     help="regenerate an Mk-model simulation instance (as simulation.py)")
    ap.add_argument("--costs", help="prediction cost CSV, same shape as --matrix (default: all 1)")
    ap.add_argument("--r01_h", type=float, default=0.5)
    ap.add_argument("--r10_h", type=float, default=0.5)
    ap.add_argument("--r01_p", type=float, default=0.5)
    ap.add_argument("--r10_p", type=float, default=0.5)
    ap.add_argument("--scale", type=float, default=100.0,
                    help="branch-length scale for the edge weights (simulations use 100)")
    ap.add_argument("--seed", type=int, default=42, help="--simulate: random seed")
    ap.add_argument("--r01", type=float, default=0.5, help="--simulate: 0->1 rate on both trees")
    ap.add_argument("--r10", type=float, default=0.5, help="--simulate: 1->0 rate on both trees")
    ap.add_argument("--corrupt", type=int, default=1000,
                    help="--simulate: number of true interactions to hide")
    ap.add_argument("--genomic_costs", action="store_true",
                    help="--simulate: genomic-informed costs (simulation.py --mode CG)")
    ap.add_argument("--elbow_grid", type=int, default=50,
                    help="pick the elbow (kneedle) from the solutions at N evenly spaced lambda "
                         "in [0,1]. The default N=50 is the sweep used for the results in the "
                         "paper and reproduces its elbows; N=0 applies kneedle to every matrix "
                         "on the front, which can select a different elbow")
    ap.add_argument("--check_grid", type=int, default=0,
                    help="verify against solve_network_cut at N evenly spaced lambda in [0,1]")
    ap.add_argument("--no_plot", action="store_true")
    ap.add_argument("--outdir", default="parametric_out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    hidden_cells = []
    if args.simulate:
        out, host_W, par_W, costs, hidden_cells = simulated_instance(args)
    else:
        matrix = pd.read_csv(args.matrix, index_col=0)
        out = build_out(args.host_tree, args.virus_tree, matrix)
        host_W, par_W = weight_matrices(out, args.r01_h, args.r10_h,
                                        args.r01_p, args.r10_p, args.scale)
        if args.costs:
            cdf = pd.read_csv(args.costs, index_col=0)
            cdf.index = cdf.index.astype(str)
            cdf.columns = cdf.columns.astype(str)
            costs = cdf.loc[out["par_leaves"], out["host_leaves"]].to_numpy(dtype=float)
        else:
            costs = build_flip_cost_matrix(len(out["par_leaves"]), len(out["host_leaves"]))

    solutions = enumerate_pareto_parametric(out, host_W, par_W, costs, verbose=True)
    for s in solutions:
        s["parsimony"] = unweighted_parsimony(out, s["new_cell_state"])
    print(f"[parametric] {len(solutions)} Pareto optimal matrices on the lower convex hull")

    front = pd.DataFrame([{
        "lambda_low": s["lambda_low"], "lambda_high": s["lambda_high"],
        "mu_low": s["mu_low"], "mu_high": s["mu_high"],
        "n_flips": len(s["flips"]), "prediction_cost": s["prediction_cost"],
        "weighted_parsimony": s["weighted_parsimony"], "parsimony": s["parsimony"],
    } for s in solutions])
    front.to_csv(os.path.join(args.outdir, "pareto_front.csv"), index=False)

    if args.elbow_grid > 0:
        lambdas = np.linspace(0.0, 1.0, args.elbow_grid)
        picked = [solution_at_lambda(solutions, lam) for lam in lambdas]
        k = find_elbow([(len(s["flips"]), s["parsimony"]) for s in picked])
        elbow = picked[k]
    else:
        k = find_elbow([(len(s["flips"]), s["parsimony"]) for s in solutions])
        elbow = solutions[k]
    print(f"[elbow] {len(elbow['flips'])} predicted interactions, "
          f"lambda in [{elbow['lambda_low']:.6g}, {elbow['lambda_high']:.6g}]")

    pd.DataFrame(elbow["new_matrix"], index=out["par_leaves"], columns=out["host_leaves"]) \
        .to_csv(os.path.join(args.outdir, "elbow_matrix.csv"))
    pd.DataFrame([(p, h) for p, h, _, _ in elbow["flips"]], columns=["virus", "host"]) \
        .to_csv(os.path.join(args.outdir, "elbow_predictions.csv"), index=False)

    summary = {
        "n_solutions": len(solutions),
        "elbow_n_flips": len(elbow["flips"]),
        "elbow_lambda_low": elbow["lambda_low"],
        "elbow_lambda_high": elbow["lambda_high"],
        "elbow_from": f"grid of {args.elbow_grid}" if args.elbow_grid > 0 else "full front",
    }
    if hidden_cells:
        from simulation_utils import compute_metrics
        summary["elbow_metrics"] = compute_metrics(hidden_cells, elbow["flips"])
        print("[elbow] metrics:", summary["elbow_metrics"])

    if args.check_grid > 0:
        # A different matrix with the same objective value is a tie: lambda lies
        # on a breakpoint, both matrices are optimal there, and which one a
        # solver returns depends on floating-point rounding.
        mismatches = ties = 0
        for lam in np.linspace(0.0, 1.0, args.check_grid):
            ref = solve_network_cut(out, host_W_matrices=host_W, par_W_matrices=par_W,
                                    flip_cost_matrix=costs, lambda_param=lam)
            ours = solution_at_lambda(solutions, lam)
            if np.array_equal(ref["new_matrix"], ours["new_matrix"]):
                status = "match"
            elif abs(ref["cut_value"] - objective(ours, lam)) <= 1e-9 * max(1.0, abs(ref["cut_value"])):
                status, ties = "tie", ties + 1
            else:
                status, mismatches = "MISMATCH", mismatches + 1
            print(f"[check] lambda={lam:.4f} flips={len(ref['flips'])} {status}")
        summary["check_grid"] = {"n": args.check_grid, "ties": ties, "mismatches": mismatches}
        print(f"[check] {args.check_grid} values of lambda: {ties} ties, {mismatches} mismatches")

    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 7))
        plt.plot(front["n_flips"], front["parsimony"], "bo-", label="Pareto front")
        plt.plot(len(elbow["flips"]), elbow["parsimony"], marker="*", markersize=20,
                 markerfacecolor="red", markeredgecolor="k", linestyle="none",
                 label="CoEvoLink solution")
        plt.xlabel("Prediction cost", fontsize=18)
        plt.ylabel("Parsimony score", fontsize=18)
        plt.legend(fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "pareto_front.svg"))
        plt.close()


if __name__ == "__main__":
    main()
