import numpy as np
from itertools import chain, combinations

def get_adj_matrix(n, edges):
    """Creates an adjacency matrix from a list of edges."""
    adj = np.zeros((n, n), dtype=int)
    for u, v in edges:
        # Adjust for 1-based indexing if necessary
        adj[u-1, v-1] += 1
    return adj

def powerset(iterable):
    """powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def solve_loopless_cycle_cover_parity(n, edges):
    """
    Calculates the parity of the number of loopless cycle covers for a given graph.
    
    This function implements the inclusion-exclusion principle.
    WARNING: This has exponential complexity and is not efficient for large graphs.
    """
    M = get_adj_matrix(n, edges)

    # Find symmetric edges for potential 2-cycles
    symmetric_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] > 0 and M[j, i] > 0:
                symmetric_edges.append(frozenset({i, j}))

    # Find all matchings in the graph of symmetric edges
    # A matching is a set of disjoint symmetric edges
    all_matchings = []
    for subset in powerset(symmetric_edges):
        is_matching = True
        nodes_in_matching = set()
        for edge in subset:
            if not nodes_in_matching.isdisjoint(edge):
                is_matching = False
                break
            nodes_in_matching.update(edge)
        if is_matching:
            all_matchings.append(subset)
    
    total_parity = 0
    print("Applying Inclusion-Exclusion Principle (mod 2):")
    # Equation: Parity = SUM_m [ product_{{i,j} in m} (M_ij * M_ji) * det(M_{V-V(m)}) ] (mod 2)
    
    for m in all_matchings:
        # Determine the sign based on the size of the matching. Modulo 2, (-1)^k is always 1.
        # So we just sum up the terms.
        
        # Calculate the weight for the matching m
        weight = 1
        for edge in m:
            u, v = tuple(edge)
            weight = (weight * M[u, v] * M[v, u])
        
        if weight == 0:
            continue

        # Get the vertices involved in the matching
        V_m = set()
        for edge in m:
            V_m.update(edge)
            
        # Create the submatrix M_{V \ V(m)}
        remaining_indices = [i for i in range(n) if i not in V_m]
        sub_M = M[np.ix_(remaining_indices, remaining_indices)]
        
        # Calculate determinant of the submatrix modulo 2
        if sub_M.shape[0] == 0:
            det_sub_M = 1
        else:
            # Using numpy to calculate determinant, then take mod 2
            det_val = np.linalg.det(sub_M)
            det_sub_M = int(round(det_val)) % 2

        term = (weight % 2 * det_sub_M) % 2
        total_parity = (total_parity + term) % 2
        
        # Outputting the calculation for each term in the sum
        v_m_str = "{" + ", ".join(map(str, sorted([v+1 for v in V_m]))) + "}"
        print(f"  Matching m = { {f'{{{u+1},{v+1}}}' for u,v in [tuple(e) for e in m]} or '{}' }: "
              f"weight = {weight}, det(M_V\\{v_m_str}) = {det_sub_M}. Term (mod 2) = {term}")


    print(f"\nFinal parity of loopless cycle covers is {total_parity}")


if __name__ == '__main__':
    # Example graph where the algorithm from thought process failed.
    # V = {1, 2, 3}, E = {(1,2), (2,1), (2,3), (3,1)}
    n_nodes = 3
    edge_list = [(1, 2), (2, 1), (2, 3), (3, 1)]
    solve_loopless_cycle_cover_parity(n_nodes, edge_list)
