import itertools

def determinant_mod2(matrix):
    """
    Computes the determinant of a matrix over GF(2) using Gaussian elimination.
    """
    n = len(matrix)
    mat = [row[:] for row in matrix]  # Make a copy

    det = 1
    for i in range(n):
        # Find pivot
        pivot_row = i
        while pivot_row < n and mat[pivot_row][i] == 0:
            pivot_row += 1

        if pivot_row == n:
            return 0  # No pivot in this column, matrix is singular

        # Swap rows
        mat[i], mat[pivot_row] = mat[pivot_row], mat[i]
        if i != pivot_row:
            det = (det + 1) % 2 # In GF(2), determinant sign doesn't matter, but for general fields it does.

        # Eliminate other rows
        for j in range(i + 1, n):
            if mat[j][i] == 1:
                for k in range(i, n):
                    mat[j][k] = (mat[j][k] + mat[i][k]) % 2
    
    # The determinant of a triangular matrix is the product of its diagonal elements.
    # After our elimination, the diagonal is all 1s if it was non-singular.
    return 1


def get_submatrix(matrix, rows_to_remove, cols_to_remove):
    """
    Creates a submatrix by removing specified rows and columns.
    """
    sub = []
    for r_idx, row in enumerate(matrix):
        if r_idx in rows_to_remove:
            continue
        new_row = []
        for c_idx, val in enumerate(row):
            if c_idx in cols_to_remove:
                continue
            new_row.append(val)
        sub.append(new_row)
    return sub

def solve_loopless_cycle_cover_parity(graph_matrix):
    """
    Calculates the parity of the number of loopless cycle covers for a graph.
    """
    n = len(graph_matrix)
    # Convert to a matrix over GF(2)
    A = [[val % 2 for val in row] for row in graph_matrix]

    # Term 1: det(A)
    det_A = determinant_mod2(A)
    
    # Term 2: Sum over all potential 2-cycles
    sum_term = 0
    
    print("Calculating parity of loopless cycle covers:")
    print(f"Number of vertices n = {n}")
    print(f"Adjacency matrix mod 2, A = {A}")
    print(f"First term, det(A) mod 2 = {det_A}")
    print("\nCalculating second term: sum_{i<j} A[i][j]*A[j][i]*det(A_ij) mod 2")
    
    for i, j in itertools.combinations(range(n), 2):
        # Check for a 2-cycle between i and j
        if A[i][j] == 1 and A[j][i] == 1:
            submatrix_ij = get_submatrix(A, {i, j}, {i, j})
            det_submatrix = determinant_mod2(submatrix_ij)
            print(f"  - Found 2-cycle between nodes {i} and {j}.")
            print(f"    Submatrix A_{{{i},{j}}} is {submatrix_ij}")
            print(f"    det(A_{{{i},{j}}}) mod 2 = {det_submatrix}")
            sum_term = (sum_term + det_submatrix) % 2
            
    print(f"\nSecond term (sum) = {sum_term}")
    
    result = (det_A + sum_term) % 2
    
    print("\nFinal Calculation:")
    print(f"Parity = (det(A) + sum) mod 2")
    print(f"Parity = ({det_A} + {sum_term}) mod 2")
    print(f"Parity = {(det_A + sum_term)} mod 2")
    print(f"Parity = {result}")

    return result

# Example graph: A 4-cycle 0->1->2->3->0 and a 2-cycle 0<->2
# This graph has one loopless cycle cover: the 4-cycle itself.
# So the parity should be 1.
# G = [[0, 1, 1, 0],
#      [0, 0, 1, 0],
#      [1, 0, 0, 1],
#      [1, 0, 0, 0]]
# Let's use a simpler example from Koiran et al.
# G = [[0, 1, 1], [1, 0, 1], [1, 1, 0]] (K3 with all edges)
# Cycle covers: (0->1->2->0) and (0->2->1->0). Both are loopless. Total 2. Parity 0.
# Another: 2-cycle (0->1->0) and self-loop on 2. Not a cycle cover of K3.
G = [[0, 1, 1, 0],
     [1, 0, 0, 1],
     [1, 0, 0, 1],
     [0, 1, 1, 0]]

# In this graph G, some loopless cycle covers are:
# (0->1->3->2->0)
# The algorithm will calculate the parity.
solve_loopless_cycle_cover_parity(G)
