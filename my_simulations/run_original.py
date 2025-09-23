import numpy as np
from Bio import Phylo

# ---------------------------
# Read NEXUS matrix with host tip names as columns
# ---------------------------
def read_nexus_with_host_columns(nexus_path, host_tip_names):
    """
    Reads a NEXUS matrix where rows = parasites, columns = hosts.
    Returns:
        parasites: list of parasite names (rows)
        hosts: list of host names (columns, from host_tip_names)
        matrix: dict[(parasite, host)] -> 0/1
    """
    matrix = {}
    parasites = []
    n_hosts = len(host_tip_names)
    
    with open(nexus_path) as f:
        lines = f.readlines()
    start_idx = 0
    for i, l in enumerate(lines):
        if l.strip().lower().startswith("matrix"):
            start_idx = i+1
            break
    
    for l in lines[start_idx:]:
        if l.strip() == ";" or l.strip().lower().startswith("end"):
            break
        parts = l.strip().split()
        if len(parts) < 2:
            continue
        parasite = parts[0]
        seq = parts[1]
        parasites.append(parasite)
        if len(seq) != n_hosts:
            raise ValueError(f"Parasite {parasite} sequence length {len(seq)} != number of hosts {n_hosts}")
        for j, c in enumerate(seq):
            h_name = host_tip_names[j]
            if c == '0':
                matrix[(parasite, h_name)] = 0
            elif c in ('1','2'):
                matrix[(parasite, h_name)] = 1
            else:  # '?' or '-'
                matrix[(parasite, h_name)] = 0
    return parasites, host_tip_names, matrix

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

# ---------------------------
# Example usage
# ---------------------------
if __name__ == "__main__":

    score = total_parsimony(
        host_tree_file="angio_25tips_origin.phy",
        parasite_tree_file="Nymphalini.phy",
        nexus_file="nymphalini.3s.nex"
    )
    print("Total parsimony (host + parasite) =", score)
