import networkx as nx
from itertools import combinations
from math import pow

def get_common_neighbors(G, node_list):
    """Computes the set of common neighbors for a list of nodes."""
    if not node_list:
        return set(G.nodes())
    
    # Start with the neighbors of the first node
    common_neighbors = set(G.neighbors(node_list[0]))
    
    # Intersect with the neighbors of the other nodes
    for i in range(1, len(node_list)):
        common_neighbors.intersection_update(G.neighbors(node_list[i]))
    return common_neighbors

def solve_count_ans_pie(G, k):
    """
    Calculates the number of answers for phi_k in G using the
    Principle of Inclusion-Exclusion. This method is exact but not efficient.
    """
    
    nodes = list(G.nodes())
    n = G.number_of_nodes()
    
    print(f"Solving for a graph with {n} nodes and k={k}.")
    print("The number of answers is given by the Inclusion-Exclusion formula:")
    print("Count = S_1 - S_2 + S_3 - ... + (-1)^(n-1) * S_n")
    print("where S_j = sum_{Y subset V, |Y|=j} |CommonNeighbors(Y)|^k")
    print("-" * 30)

    total_answers = 0
    equation_parts = []

    for j in range(1, n + 1):
        # Calculate S_j
        s_j = 0
        for y_set in combinations(nodes, j):
            # For each subset Y of size j, find common neighbors and add |CN(Y)|^k
            cn = get_common_neighbors(G, list(y_set))
            term = int(pow(len(cn), k))
            s_j += term
        
        # Apply inclusion-exclusion principle
        if j % 2 == 1:
            total_answers += s_j
            equation_parts.append(f"+ {s_j}")
        else:
            total_answers -= s_j
            equation_parts.append(f"- {s_j}")
    
    # Format the equation string for printing
    equation_str = " ".join(equation_parts)
    if equation_str.startswith("+ "):
        equation_str = equation_str[2:]
        
    print("The full calculation is:")
    print(f"{equation_str} = {total_answers}")
    print("-" * 30)
    print(f"The total number of answers is: {total_answers}")


# --- Example Usage ---
# Let's use a 3-star graph: V={0,1,2,3}, with 0 as the center.
# The common neighbor of {1,2,3} is {0}.
# The common neighbor of any other set of leaves is also {0}.
# The neighbors of the center {0} are {1,2,3}.
G_star = nx.Graph()
G_star.add_edges_from([(0, 1), (0, 2), (0, 3)])
k_val = 2

solve_count_ans_pie(G_star, k_val)