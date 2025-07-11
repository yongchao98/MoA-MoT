import numpy as np

def get_det_mod_2(matrix):
    """
    Computes the determinant of a matrix over GF(2).
    """
    mat = np.array(matrix, dtype=int)
    n = mat.shape[0]
    if n == 0:
        return 1

    # Use Gaussian elimination
    for i in range(n):
        # Find a pivot
        pivot_row = i
        while pivot_row < n and mat[pivot_row, i] == 0:
            pivot_row += 1

        if pivot_row == n:
            # No pivot in this column
            return 0

        # Swap rows to bring pivot to diagonal
        if pivot_row != i:
            mat[[i, pivot_row]] = mat[[pivot_row, i]]

        # Eliminate other 1s in this column
        for j in range(i + 1, n):
            if mat[j, i] == 1:
                mat[j, :] = (mat[j, :] + mat[i, :]) % 2

    # Determinant is the product of diagonal elements, which is always 1
    # since we ensured all pivots are 1.
    return 1

def solve_loopless_cycle_cover(M):
    """
    Calculates the parity of the number of loopless cycle covers
    for a directed multigraph given by its adjacency matrix M.

    This implementation has exponential time complexity but demonstrates
    the combinatorial approach to solving the problem.
    """
    n = len(M)
    M_mod_2 = [[M[i][j] % 2 for j in range(n)] for i in range(n)]

    # 1. Construct the conflict graph G_2
    g2_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if M[i][j] > 0 and M[j][i] > 0:
                g2_edges.append((i, j))

    # 2. Recursively generate all matchings of G_2
    # and sum determinants
    memo_det = {}
    
    # Final equation will be printed term by term
    print("The parity is calculated by summing determinants for each matching in the conflict graph.")
    print("Equation: |LCC| mod 2 = ( sum_S det(M_V(S)) ) mod 2")
    
    total_det_sum_mod_2 = 0
    
    memo_matchings = {}

    def generate_matchings_and_sum_dets(edge_index, covered_vertices):
        nonlocal total_det_sum_mod_2
        
        covered_tuple = tuple(sorted(list(covered_vertices)))
        if (edge_index, covered_tuple) in memo_matchings:
            return memo_matchings[(edge_index, covered_tuple)]

        # Base case: all edges considered for this path
        if edge_index == len(g2_edges):
            # A matching is found (defined by covered_vertices)
            # Create submatrix by removing covered vertices
            all_vertices = set(range(n))
            uncovered_vertices = sorted(list(all_vertices - covered_vertices))
            
            uncovered_tuple = tuple(uncovered_vertices)
            
            if uncovered_tuple in memo_det:
                det_val = memo_det[uncovered_tuple]
            else:
                sub_matrix = [[M_mod_2[r][c] for c in uncovered_vertices] for r in uncovered_vertices]
                det_val = get_det_mod_2(sub_matrix)
                memo_det[uncovered_tuple] = det_val
            
            print(f"  Matching covering {{{', '.join(map(str, sorted(list(covered_vertices))))}}}: det(M_V(S)) = {det_val}")
            
            return det_val

        # Recursive step
        # Option 1: Don't include the current edge in the matching
        res = generate_matchings_and_sum_dets(edge_index + 1, covered_vertices)

        # Option 2: Include the current edge if possible
        u, v = g2_edges[edge_index]
        if u not in covered_vertices and v not in covered_vertices:
            new_covered = covered_vertices.copy()
            new_covered.add(u)
            new_covered.add(v)
            res = (res + generate_matchings_and_sum_dets(edge_index + 1, new_covered)) % 2
            
        memo_matchings[(edge_index, covered_tuple)] = res
        return res
    
    total_det_sum_mod_2 = generate_matchings_and_sum_dets(0, set())

    print(f"\nFinal sum of determinants modulo 2: {total_det_sum_mod_2}")
    return total_det_sum_mod_2

if __name__ == '__main__':
    # Example Graph G=(V,E) where V={0,1,2,3} and
    # E has arcs creating two 2-cycles: (0,1),(1,0) and (2,3),(3,2)
    # The adjacency matrix M:
    n_example = 4
    M_example = [
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]
    ]

    print("--- Example ---")
    print(f"Solving for a graph with {n_example} vertices.")
    print("Adjacency Matrix M:")
    for row in M_example:
        print(row)
    
    print("\nCalculation:")
    result = solve_loopless_cycle_cover(M_example)
    print(f"\nResult: The parity of the number of loopless cycle covers is {result}.")
    print("\nAs reasoned in the explanation, the problem is NP-hard. This algorithm runs in exponential time.")

<<<A>>>