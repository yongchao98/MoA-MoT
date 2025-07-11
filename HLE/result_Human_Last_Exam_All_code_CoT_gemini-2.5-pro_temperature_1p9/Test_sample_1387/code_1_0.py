import copy

def det_gf2(matrix):
    """
    Computes the determinant of a square matrix over GF(2).
    """
    n = len(matrix)
    mat = copy.deepcopy(matrix) # Work on a copy

    det = 1
    for i in range(n):
        # Find a pivot
        pivot_row = -1
        for j in range(i, n):
            if mat[j][i] == 1:
                pivot_row = j
                break

        if pivot_row == -1:
            return 0  # No pivot in this column, determinant is 0

        # Swap rows to bring pivot to mat[i][i]
        mat[i], mat[pivot_row] = mat[pivot_row], mat[i]
        if i != pivot_row:
            det *= -1 # In GF(2), -1 is 1, so this has no effect, but good practice.
                      # We can ignore det since we only need the final parity 0 or 1.

        # Eliminate other 1s in the current column
        for j in range(n):
            if i != j and mat[j][i] == 1:
                # Add row i to row j (XOR in GF(2))
                for k in range(i, n):
                    mat[j][k] = mat[j][k] ^ mat[i][k]
    
    # After Gaussian elimination, the determinant is the product of diagonal elements.
    # In our case, the diagonal is all 1s unless the matrix was singular.
    final_det = 1
    for i in range(n):
        final_det *= mat[i][i]
        
    return final_det

def solve_loopless_cycle_cover_parity(G_adj_matrix):
    """
    Calculates the parity of the number of loopless cycle covers.
    
    Args:
    G_adj_matrix (list of lists): The adjacency matrix of the directed multigraph G.
                                 G_adj_matrix[i][j] is the number of arcs from i to j.
    """
    n = len(G_adj_matrix)
    
    # 1. Create adjacency matrix A' over GF(2)
    A_prime = [[G_adj_matrix[i][j] % 2 for j in range(n)] for i in range(n)]
    
    # 2. Create matrix C for potential 2-cycles
    C = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                C[i][j] = A_prime[i][j] * A_prime[j][i]

    # 3. Create the 2n x 2n matrix M = [[C, A'], [A'^T, 0]]
    A_prime_T = [[A_prime[j][i] for j in range(n)] for i in range(n)]
    
    M = [[0] * (2 * n) for _ in range(2 * n)]
    
    for r in range(n):
        for c in range(n):
            # Top-left block C
            M[r][c] = C[r][c]
            # Top-right block A'
            M[r][c + n] = A_prime[r][c]
            # Bottom-left block A'^T
            M[r + n][c] = A_prime_T[r][c]
            # Bottom-right block is already zero
    
    # 4. Compute the determinant of M over GF(2)
    parity = det_gf2(M)
    
    # We should output each number of the equation for clarity
    print(f"Parity of Loopless Cycle Covers of G is det(M) mod 2.")
    print("M is a 2n x 2n matrix constructed from A' (adjacency matrix mod 2) and C (2-cycle indicator matrix).")
    print(f"det(M) mod 2 = {parity}")


if __name__ == '__main__':
    # Example usage:
    # A directed 4-cycle: 0->1, 1->2, 2->3, 3->0.
    # One loopless cycle cover: the cycle itself. Parity = 1.
    print("Example 1: A directed 4-cycle")
    adj_matrix1 = [
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0]
    ]
    solve_loopless_cycle_cover_parity(adj_matrix1)
    print("-" * 20)

    # Two disjoint 2-cycles: 0<->1, 2<->3
    # One cycle cover exists, but it's not loopless. Parity = 0.
    print("Example 2: Two disjoint 2-cycles")
    adj_matrix2 = [
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]
    ]
    solve_loopless_cycle_cover_parity(adj_matrix2)
    print("-" * 20)
    
    # A 3-cycle 0->1->2->0 with a reversing edge 1->0
    # The only cycle cover is the 3-cycle itself, which is loopless. Parity = 1.
    print("Example 3: A 3-cycle with a reversing edge")
    adj_matrix3 = [
        [0, 1, 0],
        [1, 0, 1],
        [0, 0, 0]
    ]
    # To make it a cycle cover, vertex 2 needs an outgoing edge and vertex 0 an incoming one
    # let's use 2->0.
    adj_matrix3[2][0] = 1
    solve_loopless_cycle_cover_parity(adj_matrix3)
