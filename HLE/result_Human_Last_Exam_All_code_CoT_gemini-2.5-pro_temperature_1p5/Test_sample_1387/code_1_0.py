import numpy as np
from itertools import combinations

def get_submatrix(matrix, indices_to_remove):
    """
    Removes specified rows and columns from a matrix.
    """
    if not indices_to_remove:
        return matrix
    # Create a mask for rows/columns to keep
    mask = np.ones(matrix.shape[0], dtype=bool)
    mask[list(indices_to_remove)] = False
    return matrix[mask, :][:, mask]

def calculate_determinant_mod2(matrix):
    """
    Calculates the determinant of a matrix modulo 2.
    """
    if matrix.size == 0:
        return 1
    det = np.linalg.det(matrix)
    return int(round(det)) % 2

def generate_matchings(edges):
    """
    Recursively generates all matchings of a graph given its edge list.
    """
    if not edges:
        yield []
        return

    first_edge, *rest_edges = edges
    u, v = first_edge

    # Case 1: The first edge is not in the matching
    for matching in generate_matchings(rest_edges):
        yield matching

    # Case 2: The first edge is in the matching
    # Remove all other edges incident to u or v
    compatible_edges = [
        (n1, n2) for n1, n2 in rest_edges if n1 != u and n1 != v and n2 != u and n2 != v
    ]
    for matching in generate_matchings(compatible_edges):
        yield [first_edge] + matching

def solve_loopless_cycle_cover_parity(adj_matrix):
    """
    Calculates the parity of the number of loopless cycle covers
    using the inclusion-exclusion principle. This implementation is
    correct but not polynomial-time in the worst case.
    
    Args:
        adj_matrix (np.ndarray): The adjacency matrix of the graph.
    """
    n = adj_matrix.shape[0]
    adj_matrix_mod2 = adj_matrix % 2
    
    # Find all potential 2-cycles
    two_cycle_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if adj_matrix_mod2[i, j] == 1 and adj_matrix_mod2[j, i] == 1:
                two_cycle_edges.append((i, j))
    
    print(f"Graph has {n} vertices.")
    if two_cycle_edges:
        print(f"Found potential 2-cycles between vertex pairs: {two_cycle_edges}")
    else:
        print("No potential 2-cycles found.")
    
    total_parity = 0
    
    print("\nCalculating terms using Inclusion-Exclusion (sum over all 2-cycle matchings):")
    
    for matching in generate_matchings(two_cycle_edges):
        vertices_in_matching = set()
        for u, v in matching:
            vertices_in_matching.add(u)
            vertices_in_matching.add(v)
            
        sub_matrix = get_submatrix(adj_matrix_mod2, list(vertices_in_matching))
        det_mod2 = calculate_determinant_mod2(sub_matrix)
        
        term_str = f"det(A_V"
        if vertices_in_matching:
            term_str += f"\\{{{','.join(map(str, sorted(list(vertices_in_matching))))}\\}}}"
        term_str += f") mod 2 = {det_mod2}"
        print(term_str)
        
        total_parity = (total_parity + det_mod2) % 2

    print(f"\nFinal sum modulo 2 is: {total_parity}")
    return total_parity

if __name__ == '__main__':
    # Example Graph: A directed graph with 4 vertices.
    # It consists of two disjoint 2-cycles: 0 <-> 1 and 2 <-> 3.
    # It has one cycle cover: {(0,1), (1,0), (2,3), (3,2)}.
    # This cover is NOT loopless.
    # So, the number of loopless cycle covers is 0, and the parity is 0.
    
    G_adj = np.array([
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]
    ])
    
    solve_loopless_cycle_cover_parity(G_adj)
