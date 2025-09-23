import math, random
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import copy
import networkx as nx
from matplotlib.patches import Patch
from Bio.Nexus.Nexus import Nexus
from kneed import KneeLocator
from scipy.signal import savgol_filter
import os

# ---------------------------
# Compute log-likelihood of a simulation using assigned internal states
# ---------------------------
def compute_current_ll(tree, leaf_states, internal_states, r01=0.5, r10=0.5, eps=1e-300):
    """
    Compute log-likelihood along the tree using exact assigned states
    internal_states: dict {clade: state} including leaves
    """
    loglik = 0.0
    for parent in tree.find_clades(order="preorder"):
        parent_state = internal_states[parent]
        for child in parent.clades:
            child_state = internal_states[child]
            t = child.branch_length if child.branch_length else 0.0
            P = expQt_2state(t, r01, r10)
            prob = P[parent_state, child_state]
            loglik += math.log(prob + eps)
    return loglik

# ---------------------------
# Total log-likelihood using assigned states
# ---------------------------
def total_current_ll(host_trees, par_trees, cell_state, host_leaves, par_leaves, r01=0.5, r10=0.5):
    total_ll = 0.0
    # Host tree: each parasite row = character
    for i, host_tree in enumerate(host_trees):
        for p in par_leaves:
            leaf_states = {h: cell_state[(p, h)] for h in host_leaves}
            # Internal states = the states assigned during simulation
            # host_tree already contains internal nodes assigned in simulation
            internal_states = host_tree.__dict__.get("states", {})  # optional: store during simulation
            if not internal_states:
                # fallback: only leaf states
                internal_states = {h_node: leaf_states[h_node.name] for h_node in host_tree.get_terminals()}
            total_ll += compute_current_ll(host_tree, leaf_states, internal_states, r01, r10)
    # Parasite tree: each host column = character
    for j, par_tree in enumerate(par_trees):
        for h in host_leaves:
            leaf_states = {p: cell_state[(p, h)] for p in par_leaves}
            internal_states = par_tree.__dict__.get("states", {})
            if not internal_states:
                internal_states = {p_node: leaf_states[p_node.name] for p_node in par_tree.get_terminals()}
            total_ll += compute_current_ll(par_tree, leaf_states, internal_states, r01, r10)
    return total_ll


# ---------------------------
# Sankoff parsimony
# ---------------------------
def sankoff(tree, leaf_states, cost_matrix):
    k = cost_matrix.shape[0]
    def recurse(clade):
        if clade.is_terminal():
            state = leaf_states[clade.name]
            cost = np.full(k, np.inf)
            cost[state] = 0
            return cost
        child_costs = [recurse(c) for c in clade.clades]
        total_cost = np.zeros(k)
        for s in range(k):
            total_cost[s] = sum(
                min(child[s2] + cost_matrix[s, s2] for s2 in range(k))
                for child in child_costs
            )
        return total_cost
    root_costs = recurse(tree.root)
    return min(root_costs)

# ---------------------------
# Total parsimony (host + parasite)
# ---------------------------
def total_parsimony(host_tree_file, parasite_tree_file, nexus_file=None, cell_state=None):
    host_tree = Phylo.read(host_tree_file, "newick")
    parasite_tree = Phylo.read(parasite_tree_file, "newick")
    host_tip_names = [t.name for t in host_tree.get_terminals()]
    parasite_tip_names = [t.name for t in parasite_tree.get_terminals()]
    
    if cell_state is None and nexus_file is not None:
        parasites, hosts, cell_state = read_nexus_with_host_columns(nexus_file, host_tip_names)
    elif cell_state is None:
        raise ValueError("Either nexus_file or cell_state must be provided")
    else:
        parasites = [p for (p,h) in cell_state.keys() if h == host_tip_names[0]]
    
    C = np.array([[0,1],[1,0]])
    total_score = 0
    
    # Host tree: each parasite row = character
    for p in parasites:
        leaf_states = {h: cell_state[(p,h)] for h in host_tip_names}
        total_score += sankoff(host_tree, leaf_states, C)
    
    # Parasite tree: each host column = character
    for h in host_tip_names:
        leaf_states = {p: cell_state[(p,h)] for p in parasites}
        total_score += sankoff(parasite_tree, leaf_states, C)
    
    return total_score



def avg_P_same(t, r01=0.5, r10=0.5):
    r_sum = r01 + r10
    if r_sum == 0.0:
        return 1.0
    e = math.exp(-r_sum * t)
    num = (r10*r10 + r01*r01) + 2.0 * r01 * r10 * e
    return num / (r_sum * r_sum)

# ---------------------------
# Mk(2) transition matrix
# ---------------------------
def expQt_2state(t, r01=0.5, r10=0.5):
    """Transition probability matrix for 2 states with asymmetric rates"""
    r_sum = r01 + r10
    e = math.exp(-r_sum * t)
    
    p00 = r10/r_sum + r01/r_sum * e
    p01 = r01/r_sum - r01/r_sum * e
    p10 = r10/r_sum - r10/r_sum * e
    p11 = r01/r_sum + r10/r_sum * e
    
    return np.array([[p00, p01],
                     [p10, p11]])
  # rows: from-state 0/1, cols: to-state 0/1


def build_weight_matrix(tree, r01, r10):
    """
    Create a weight matrix for a tree.
    All edges use the same transition probability (avg_P_same).
    """
    
    # Create adjacency matrix
    leaves = tree.get_terminals()
    nodes = list(tree.find_clades(order="level"))
    n = len(nodes)
    W = np.zeros((n, n))

    for parent in nodes:
        for child in parent.clades:
            i, j = nodes.index(parent), nodes.index(child)
            W[i, j] = avg_P_same(child.branch_length or 0.0, r01, r10)
            # W[i, j] = child.branch_length
            W[j, i] = W[i, j]  # symmetric

    return W, nodes

def build_flip_cost_matrix(n_parasites, n_hosts, cost=1.0):
    """
    Create a flip cost matrix for parasite × host interactions.
    Default = all 1s.
    """
    return np.full((n_parasites, n_hosts), cost)



# ---------------------------
# Simulate states along a single tree
# ---------------------------
def simulate_internal_states(tree, r01=0.5, r10=0.5, seed=None):
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    root = tree.root
    states = {}
    states[root] = random.choice([0, 1])
    
    def recurse(parent):
        parent_state = states[parent]
        for child in parent.clades:
            t = child.branch_length if child.branch_length else 0.0
            P = expQt_2state(t, r01=r01, r10=r10)
            probs = P[parent_state] / P[parent_state].sum()
            states[child] = int(np.random.choice([0,1], p=probs))
            recurse(child)
    recurse(root)
    # store in tree object
    tree.__dict__["states"] = states
    return states

# ---------------------------
# Build parent map
# ---------------------------
def build_parent_map(tree):
    parent = {}
    for p in tree.find_clades(order="preorder"):
        for c in p.clades:
            parent[c] = p
    return parent

# ---------------------------
# Tree depth and rescaling
# ---------------------------
def get_tree_depth(tree):
    depths = tree.depths()
    return max(depths.values())

def rescale_tree(tree):
    max_depth = get_tree_depth(tree)
    if max_depth == 0:
        return tree
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length /= max_depth
    return tree

# ---------------------------
# Simulate network of m+n trees
# ---------------------------
def simulate_network(host_tree_file, parasite_tree_file, r01_p=0.5, r10_p=0.5, r01_h=0.5, r10_h=0.5, seed=None):
    base_host_tree = Phylo.read(host_tree_file, "newick")
    base_par_tree = Phylo.read(parasite_tree_file, "newick")
    
    base_host_tree = rescale_tree(base_host_tree)
    base_par_tree = rescale_tree(base_par_tree)
    
    host_leaves = list(base_host_tree.get_terminals())
    par_leaves = list(base_par_tree.get_terminals())

    m = len(host_leaves)
    n = len(par_leaves)
    
    host_trees = []
    host_states_list = []
    host_parent_list = []
    for i in range(n):
        t = copy.deepcopy(base_host_tree)
        states = simulate_internal_states(t, r01_h, r10_h, seed=seed)
        host_trees.append(t)
        host_states_list.append(states)
        host_parent_list.append(build_parent_map(t))

    par_trees = []
    par_states_list = []
    par_parent_list = []
    for i in range(m):
        t = copy.deepcopy(base_par_tree)
        states = simulate_internal_states(t, r01_p, r10_p, seed=seed)
        par_trees.append(t)
        par_states_list.append(states)
        par_parent_list.append(build_parent_map(t))
    
    # ---------------------------
    # Parasites x Hosts
    # ---------------------------
    cell_state = {}
    cell_matrix = {}
    row_idx = {(0,0):0, (0,1):1, (1,0):2, (1,1):3}

    for i, p_leaf in enumerate(par_leaves):
        for j, h_leaf in enumerate(host_leaves):
            # Find corresponding clades
            h_clade_copy = [l for l in host_trees[i].get_terminals() if l.name == h_leaf.name][0]
            p_clade_copy = [l for l in par_trees[j].get_terminals() if l.name == p_leaf.name][0]

            h_parent = host_parent_list[i][h_clade_copy]
            p_parent = par_parent_list[j][p_clade_copy]

            Ph = expQt_2state(h_clade_copy.branch_length or 0.0, r01_h, r10_h)
            Pp = expQt_2state(p_clade_copy.branch_length or 0.0, r01_p, r10_p)

            M = np.zeros((4,2))
            for a in (0,1):
                for b in (0,1):
                    r = row_idx[(a,b)]
                    M[r,0] = Ph[a,0]*Pp[b,0]
                    M[r,1] = Ph[a,1]*Pp[b,1]
                    if M[r].sum() > 0:
                        M[r] /= M[r].sum()

            h_state = host_states_list[i][h_parent]
            p_state = par_states_list[j][p_parent]
            r_star = row_idx[(h_state, p_state)]
            sampled = int(np.random.choice([0,1], p=M[r_star]))

            cell_state[(p_leaf.name, h_leaf.name)] = sampled
            cell_matrix[(p_leaf.name, h_leaf.name)] = M

    return {
        "host_trees": host_trees,
        "par_trees": par_trees,
        "cell_state": cell_state,
        "cell_matrix": cell_matrix,
        "host_leaves": [n.name for n in host_leaves],
        "par_leaves": [n.name for n in par_leaves]
    }

# ---------------------------
# Get parasite x host matrix
# ---------------------------
def get_interaction_matrix(out):
    parasites = out["par_leaves"]
    hosts = out["host_leaves"]
    mat = np.zeros((len(parasites), len(hosts)), dtype=int)
    for i, p in enumerate(parasites):
        for j, h in enumerate(hosts):
            mat[i,j] = out["cell_state"].get((p,h), 0)  # key is (parasite, host)
    return mat, parasites, hosts

# ---------------------------
# Plot matrix
# ---------------------------
# def plot_matrix(mat, row_names, col_names, filename="matrix.png", xlabel="Hosts", ylabel="Parasites"):
#     fig, ax = plt.subplots(figsize=(len(col_names)*0.6, len(row_names)*0.6))
#     cmap = mcolors.ListedColormap(['blue','red'])
#     ax.imshow(mat, cmap=cmap, vmin=0, vmax=1)
#     ax.set_xticks(np.arange(len(col_names)))
#     ax.set_yticks(np.arange(len(row_names)))
#     ax.set_xticklabels(col_names, rotation=90, fontsize=8)
#     ax.set_yticklabels(row_names, fontsize=8)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.set_title(f"{ylabel} x {xlabel} Matrix")
    
#     ax.set_xticks(np.arange(-0.5, len(col_names), 1), minor=True)
#     ax.set_yticks(np.arange(-0.5, len(row_names), 1), minor=True)
#     ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.5)
#     ax.tick_params(which="minor", bottom=False, left=False)
#     plt.tight_layout()
#     plt.savefig(filename, dpi=300)
#     plt.close(fig)
#     print(f"Matrix saved as {filename}")


def plot_matrix(mat, row_names, col_names, filename="matrix.png",
                xlabel="Hosts", ylabel="Parasites",
                highlight=None):
    """
    Plot binary matrix with optional highlighted cells.
    
    highlight: dict with keys as category name and values as 
               list of (row_name, col_name) coordinates.
               Example:
               {
                   "corrupted": [(p1,h1), (p2,h3)],
                   "flipped": [(p4,h2)]
               }
    """
    fig, ax = plt.subplots(figsize=(len(col_names)*0.6, len(row_names)*0.6))

    # Base: 0=blue, 1=red
    color_mat = np.copy(mat).astype(float)

    # Reserve higher codes for highlights
    code_map = {"corrupted": 2, "flipped": 3}
    if highlight is not None:
        for category, cells in highlight.items():
            if category not in code_map:
                continue
            for p, h in cells:
                i = row_names.index(p)
                j = col_names.index(h)
                color_mat[i, j] = code_map[category]

    # Custom colormap: 0=blue, 1=red, 2=orange (corrupted), 3=green (flipped)
    cmap = mcolors.ListedColormap(['blue', 'red', 'orange', 'green'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    im = ax.imshow(color_mat, cmap=cmap, norm=norm)

    ax.set_xticks(np.arange(len(col_names)))
    ax.set_yticks(np.arange(len(row_names)))
    ax.set_xticklabels(col_names, rotation=90, fontsize=8)
    ax.set_yticklabels(row_names, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_title(f"{ylabel} x {xlabel} Matrix")

    # Grid lines
    ax.set_xticks(np.arange(-0.5, len(col_names), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_names), 1), minor=True)
    ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # legend_elements = [
    #     Patch(facecolor='blue', label='0'),
    #     Patch(facecolor='red', label='1'),
    #     Patch(facecolor='orange', label='Randomly corrupted (1→0)'),
    #     Patch(facecolor='green', label='Flipped by algorithm'),
    # ]
    # ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close(fig)
    print(f"Matrix saved as {filename}")

def corrupt_matrix(mat, num_flips, seed=0):
    """Randomly flip num_flips of 1→0 and return corrupted matrix + flipped indices."""
    rng = random.Random(seed)
    mat_corrupt = mat.copy()
    ones = list(zip(*np.where(mat == 1)))
    chosen = rng.sample(ones, min(num_flips, len(ones)))
    for i, j in chosen:
        mat_corrupt[i, j] = 0
    return mat_corrupt, chosen  # chosen = list of (i,j) where true 1 was hidden

def compare_predictions(true_mat, corrupt_mat, predicted_mat):
    """Classify each cell for visualization."""
    n, m = true_mat.shape
    color_mat = np.zeros((n, m), dtype=int)

    for i in range(n):
        for j in range(m):
            true_val = true_mat[i, j]
            corrupt_val = corrupt_mat[i, j]
            pred_val = predicted_mat[i, j]

            if true_val == 1 and corrupt_val == 0 and pred_val == 1:
                color_mat[i, j] = 2  # correctly recovered (green)
            elif pred_val != true_val:
                color_mat[i, j] = 3  # wrong flip (orange)
            elif true_val == 0:
                color_mat[i, j] = 0  # true 0 (blue)
            else:
                color_mat[i, j] = 1  # true 1 (red)

    return color_mat

# ---------------------------
# Total log-likelihood using assigned states
# ---------------------------
def total_current_ll(host_trees, par_trees, cell_state, host_leaves, par_leaves, r01_h=0.5, r10_h=0.5, r01_p=0.5, r10_p=0.5):
    total_ll = 0.0
    # Host tree: each parasite row = character
    for i, host_tree in enumerate(host_trees):
        for p in par_leaves:
            leaf_states = {h: cell_state[(p, h)] for h in host_leaves}
            # Internal states = the states assigned during simulation
            # host_tree already contains internal nodes assigned in simulation
            internal_states = host_tree.__dict__.get("states", {})  # optional: store during simulation
            if not internal_states:
                # fallback: only leaf states
                internal_states = {h_node: leaf_states[h_node.name] for h_node in host_tree.get_terminals()}
            total_ll += compute_current_ll(host_tree, leaf_states, internal_states, r01_h, r10_h)
    # Parasite tree: each host column = character
    for j, par_tree in enumerate(par_trees):
        for h in host_leaves:
            leaf_states = {p: cell_state[(p, h)] for p in par_leaves}
            internal_states = par_tree.__dict__.get("states", {})
            if not internal_states:
                internal_states = {p_node: leaf_states[p_node.name] for p_node in par_tree.get_terminals()}
            total_ll += compute_current_ll(par_tree, leaf_states, internal_states, r01_p, r10_p)
    return total_ll


# # ---------------------------
# # Sankoff parsimony
# # ---------------------------
# def sankoff(tree, leaf_states, cost_matrix):
#     k = cost_matrix.shape[0]
#     def recurse(clade):
#         if clade.is_terminal():
#             state = leaf_states[clade.name]
#             cost = np.full(k, np.inf)
#             cost[state] = 0
#             return cost
#         child_costs = [recurse(c) for c in clade.clades]
#         total_cost = np.zeros(k)
#         for s in range(k):
#             total_cost[s] = sum(
#                 min(child[s2] + cost_matrix[s, s2] for s2 in range(k))
#                 for child in child_costs
#             )
#         return total_cost
#     root_costs = recurse(tree.root)
#     return min(root_costs)

# ---------------------------
# Total parsimony (host + parasite)
# ---------------------------
# def total_parsimony(host_tree_file, parasite_tree_file, nexus_file=None, cell_state=None):
#     host_tree = Phylo.read(host_tree_file, "newick")
#     parasite_tree = Phylo.read(parasite_tree_file, "newick")
#     host_tip_names = [t.name for t in host_tree.get_terminals()]
#     parasite_tip_names = [t.name for t in parasite_tree.get_terminals()]
    
#     if cell_state is None and nexus_file is not None:
#         parasites, hosts, cell_state = read_nexus_with_host_columns(nexus_file, host_tip_names)
#     elif cell_state is None:
#         raise ValueError("Either nexus_file or cell_state must be provided")
#     else:
#         parasites = [p for (p,h) in cell_state.keys() if h == host_tip_names[0]]
    
#     C = np.array([[0,1],[1,0]])
#     total_score = 0
    
#     # Host tree: each parasite row = character
#     for p in parasites:
#         leaf_states = {h: cell_state[(p,h)] for h in host_tip_names}
#         total_score += sankoff(host_tree, leaf_states, C)
    
#     # Parasite tree: each host column = character
#     for h in host_tip_names:
#         leaf_states = {p: cell_state[(p,h)] for p in parasites}
#         total_score += sankoff(parasite_tree, leaf_states, C)
    
#     return total_score

# ---------------------------
# Sankoff parsimony
# ---------------------------
def sankoff(tree, leaf_states, cost_matrix):
    k = cost_matrix.shape[0]
    def recurse(clade):
        if clade.is_terminal():
            state = leaf_states[clade.name]
            cost = np.full(k, np.inf)
            cost[state] = 0
            return cost
        child_costs = [recurse(c) for c in clade.clades]
        total_cost = np.zeros(k)
        for s in range(k):
            total_cost[s] = sum(
                min(child[s2] + cost_matrix[s, s2] for s2 in range(k))
                for child in child_costs
            )
        return total_cost
    root_costs = recurse(tree.root)
    return min(root_costs)

def sankoff_ml_loglik(tree, leaf_states, r01=0.5, r10=0.5, eps=1e-300):
    k = 2
    def recurse(clade):
        if clade.is_terminal():
            state = leaf_states[clade.name]
            cost = np.full(k, np.inf)
            cost[state] = 0.0
            return cost
        child_costs = [recurse(c) for c in clade.clades]
        total_cost = np.zeros(k)
        for s in range(k):
            total_cost[s] = sum(
                min(child[s2] - math.log(expQt_2state(cl.branch_length or 0, r01,r10)[s,s2]+eps)
                    for s2 in range(k))
                for cl, child in zip(clade.clades, child_costs)
            )
        return total_cost
    root_costs = recurse(tree.root)
    root_total = root_costs - math.log(0.5+eps)
    return -np.min(root_total)


# ---------------------------
# Total parsimony (host + parasite)
# ---------------------------
def total_parsimony(host_tree_file, parasite_tree_file, nexus_file=None, cell_state=None):
    host_tree = Phylo.read(host_tree_file, "newick")
    parasite_tree = Phylo.read(parasite_tree_file, "newick")
    host_tip_names = [t.name for t in host_tree.get_terminals()]
    parasite_tip_names = [t.name for t in parasite_tree.get_terminals()]
    
    if cell_state is None and nexus_file is not None:
        parasites, hosts, cell_state = read_nexus_with_host_columns(nexus_file, host_tip_names)
    elif cell_state is None:
        raise ValueError("Either nexus_file or cell_state must be provided")
    else:
        parasites = [p for (p,h) in cell_state.keys() if h == host_tip_names[0]]
    
    C = np.array([[0,1],[1,0]])
    total_score = 0
    
    # Host tree: each parasite row = character
    for p in parasites:
        leaf_states = {h: cell_state[(p,h)] for h in host_tip_names}
        total_score += sankoff(host_tree, leaf_states, C)
    
    # Parasite tree: each host column = character
    for h in host_tip_names:
        leaf_states = {p: cell_state[(p,h)] for p in parasites}
        total_score += sankoff(parasite_tree, leaf_states, C)
    
    return total_score


def solve_network_cut(
    out,
    host_W_matrices,
    par_W_matrices,
    flip_cost_matrix,
    lambda_param=0.5
):
    """
    Solve network cut using precomputed weight and flip cost matrices.

    host_W_matrices: list of adjacency-like weight matrices, one per parasite tree
    par_W_matrices: list of adjacency-like weight matrices, one per host tree
    flip_cost_matrix: 2D array (parasites × hosts), flip cost per cell
    """
    G = nx.DiGraph()
    source, sink = "SOURCE", "SINK"

    cell_state = out["cell_state"]
    hosts = out["host_leaves"]
    parasites = out["par_leaves"]

    # ---------------------------
    # Host trees: connect internal → internal, or internal → CELL
    # ---------------------------
    for p_idx, (p, t) in enumerate(zip(parasites, out["host_trees"])):
        W, nodes = host_W_matrices[p_idx]  # weight matrix + node order
        for parent in t.find_clades(order="preorder"):
            for child in parent.clades:
                w = lambda_param * W[nodes.index(parent), nodes.index(child)]
                # w = lambda_param * avg_P_same(child.branch_length or 0.0, r01=r01_h, r10=r10_h)
                if child.is_terminal():
                    node = f"CELL_{p}_{child.name}"
                    G.add_edge(parent, node, capacity=w)
                    G.add_edge(node, parent, capacity=w)
                else:
                    G.add_edge(parent, child, capacity=w)
                    G.add_edge(child, parent, capacity=w)

    # ---------------------------
    # Parasite trees: connect internal → internal, or internal → CELL
    # ---------------------------
    for h_idx, (h, t) in enumerate(zip(hosts, out["par_trees"])):
        W, nodes = par_W_matrices[h_idx]
        for parent in t.find_clades(order="preorder"):
            for child in parent.clades:
                w = lambda_param * W[nodes.index(parent), nodes.index(child)]
                # w = lambda_param * avg_P_same(child.branch_length or 0.0, r01=r01_p, r10=r10_p)
                if child.is_terminal():
                    node = f"CELL_{child.name}_{h}"
                    G.add_edge(parent, node, capacity=w)
                    G.add_edge(node, parent, capacity=w)
                else:
                    G.add_edge(parent, child, capacity=w)
                    G.add_edge(child, parent, capacity=w)

    # ---------------------------
    # Source/sink for cells (use flip_cost_matrix)
    # ---------------------------
    for p_idx, p in enumerate(parasites):
        for h_idx, h in enumerate(hosts):
            node = f"CELL_{p}_{h}"
            state = cell_state[(p, h)]
            if state == 0:
                G.add_edge(source, node,
                           capacity=(1 - lambda_param) * flip_cost_matrix[p_idx, h_idx])
            else:  # state == 1
                G.add_edge(node, sink, capacity=float("inf"))

    # ---------------------------
    # Min-cut
    # ---------------------------
    cut_value, partition = nx.minimum_cut(G, source, sink)
    reachable, non_reachable = partition

    new_cell_state = {}
    flips = []
    for (p, h), old in cell_state.items():
        node = f"CELL_{p}_{h}"
        new = 0 if node in reachable else 1
        new_cell_state[(p, h)] = new
        if new != old:
            flips.append((p, h, old, new))

    # Build new matrix
    new_matrix = np.zeros((len(parasites), len(hosts)), dtype=int)
    for i, p in enumerate(parasites):
        for j, h in enumerate(hosts):
            new_matrix[i, j] = new_cell_state[(p, h)]

    return {
        "cut_value": cut_value,
        "flips": flips,
        "new_matrix": new_matrix,
        "new_cell_state": new_cell_state
    }



def binary_search_lambda(out, parasites, hosts, hidden_cells,
                         target_flips, host_W_matrices, par_W_matrices, flip_cost_matrix, tol=3, max_iter=20):
    """
    Binary search for lambda that yields ~target_flips (within tolerance).
    Returns best lambda and the cut_result for that lambda.
    """
    low, high = 0.0, 1.0
    best_lambda, best_result = None, None    

    for it in range(max_iter):
        mid = (low + high) / 2
        cut_result = solve_network_cut(
            out, host_W_matrices=host_W_matrices, par_W_matrices=par_W_matrices,
            flip_cost_matrix=flip_cost_matrix,
            lambda_param=mid
        )
        flips = cut_result["flips"]
        num_flips = len(flips)

        print(f"[iter {it}] λ={mid:.4f} → flips={num_flips}")

        # Check tolerance condition
        if abs(num_flips - target_flips) <= tol:
            best_lambda, best_result = mid, cut_result
            break

        # Binary search update
        if num_flips > target_flips + tol:
            high = mid
        elif num_flips < target_flips - tol:
            low = mid
        else:
            best_lambda, best_result = mid, cut_result
            break

        # Keep track of closest so far
        if best_result is None or abs(num_flips - target_flips) < abs(len(best_result["flips"]) - target_flips):
            best_lambda, best_result = mid, cut_result

    # Evaluate correctness
    flips = best_result["flips"]
    flipped_cells = [(p, h) for p, h, _, _ in flips]
    correctly_recovered = sum(cell in flipped_cells for cell in hidden_cells)

    print(f"\nBest λ={best_lambda:.4f} "
          f"→ flips={len(flips)}, correct={correctly_recovered}/{target_flips}")

    # # Save visualization
    # highlight = {
    #     "corrupted": hidden_cells,
    #     "flipped": flipped_cells
    # }
    # plot_matrix(best_result["new_matrix"], parasites, hosts,
    #             filename=f"experiments/recovered_lambda{best_lambda:.4f}.png",
    #             highlight=highlight)

    return best_lambda, best_result




    



def compute_metrics(hidden_cells, flips):
    flipped_cells = [(p, h) for p, h, _, _ in flips]
    TP = sum(cell in flipped_cells for cell in hidden_cells)
    FP = sum(cell not in hidden_cells for cell in flipped_cells)
    FN = sum(cell not in flipped_cells for cell in hidden_cells)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {"precision": precision, "recall": recall, "f1": f1, "TP": TP, "FP": FP, "FN": FN}

