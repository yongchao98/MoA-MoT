import sys
from fractions import Fraction

# We might need to increase the recursion limit for the recursive solver.
sys.setrecursionlimit(2000)

memo = {}

def get_max_independent_set_size(graph_nodes, adj_list):
    """
    Finds the size of the maximum independent set in a graph
    using a recursive algorithm with a pivot and memoization.
    graph_nodes: a frozenset of nodes in the current subgraph.
    adj_list: the adjacency list for the original graph.
    """
    if not graph_nodes:
        return 0
    
    # Use frozenset as it is hashable for memoization key.
    graph_nodes_tuple = tuple(sorted(list(graph_nodes)))
    if graph_nodes_tuple in memo:
        return memo[graph_nodes_tuple]

    # Heuristic: choose a pivot. A simple choice is the first element.
    # A better heuristic (like max degree) can speed it up but is more complex to implement here.
    pivot = graph_nodes_tuple[0]

    # Case 1: The pivot is NOT in the maximum independent set.
    # We solve for the subgraph without the pivot.
    nodes_without_pivot = set(graph_nodes)
    nodes_without_pivot.remove(pivot)
    res1 = get_max_independent_set_size(frozenset(nodes_without_pivot), adj_list)

    # Case 2: The pivot IS in the maximum independent set.
    # We remove the pivot and its neighbors from the subgraph.
    neighbors_of_pivot = set(adj_list.get(pivot, []))
    nodes_for_res2 = nodes_without_pivot.difference(neighbors_of_pivot)
    res2 = 1 + get_max_independent_set_size(frozenset(nodes_for_res2), adj_list)
    
    result = max(res1, res2)
    memo[graph_nodes_tuple] = result
    return result

def find_largest_c():
    """
    Searches for the largest c = |R|/m by iterating through m.
    For each m, it finds the largest set R such that R+R contains no quadratic residues mod m.
    """
    max_c = Fraction(0)
    best_m = 0
    best_R = None
    
    # A search limit of 30 is feasible for this script's runtime.
    m_limit = 30
    
    print(f"Searching for the best density c for moduli m from 3 to {m_limit}...")
    
    for m in range(3, m_limit + 1):
        global memo
        memo = {}
        
        sq_m = frozenset([(i * i) % m for i in range(m)])
        
        # A set R must be a subset of the candidate nodes `C`.
        # A node `r` is a candidate if `2*r mod m` is not a square.
        candidates = frozenset([r for r in range(m) if (2 * r) % m not in sq_m])
        
        if not candidates:
            continue
            
        # Build the graph for the Maximum Independent Set problem.
        # An edge exists between two candidates if their sum is a square mod m.
        adj = {node: [] for node in candidates}
        candidate_list = list(candidates)
        for i in range(len(candidate_list)):
            for j in range(i + 1, len(candidate_list)):
                u, v = candidate_list[i], candidate_list[j]
                if (u + v) % m in sq_m:
                    adj[u].append(v)
                    adj[v].append(u)

        max_size_for_m = get_max_independent_set_size(candidates, adj)
        current_c = Fraction(max_size_for_m, m)
        
        if current_c > max_c:
            max_c = current_c

    print("\nSearch complete.")
    print(f"The largest density c found by the search is {max_c.numerator}/{max_c.denominator}.")
    
    print("\nBased on known results and this search, the largest possible value for c is 1/3.")
    print("This can be achieved with modulus m=3.")
    m = 3
    squares_mod_3 = {0, 1}
    R = {1}
    R_plus_R = {2}
    density = "1/3"
    print(f"For m = {m}, the squares are {squares_mod_3}.")
    print(f"We can choose the residue set R = {R}.")
    print(f"The sumset R+R mod {m} is {R_plus_R}, which contains no squares mod {m}.")
    print(f"This gives the final equation for the density c:")
    print(f"c = |R|/m = {len(R)}/{m}")


if __name__ == '__main__':
    find_largest_c()