import sys

def det_mod2(matrix):
    """
    Computes the determinant of a square matrix over F_2 (modulo 2)
    using Gaussian elimination.
    """
    n = len(matrix)
    mat = [row[:] for row in matrix] # Make a copy to modify

    det = 1
    for j in range(n): # Iterate through columns (pivots)
        # Find a row with a 1 in the current column
        pivot_row = -1
        for i in range(j, n):
            if mat[i][j] == 1:
                pivot_row = i
                break

        if pivot_row == -1:
            # If no pivot found, the column is all zeros.
            # The matrix is singular, determinant is 0.
            return 0
        
        # Swap rows to bring the pivot to the diagonal
        if pivot_row != j:
            mat[j], mat[pivot_row] = mat[pivot_row], mat[j]
            # In F_2, row swaps don't change the sign of the determinant (-1 = 1)

        # Update determinant (product of diagonal elements)
        # This is implicitly handled since we ensure diagonal is 1
        
        # Eliminate other 1s in the current column
        for i in range(n):
            if i != j and mat[i][j] == 1:
                # Add the pivot row to the current row (XOR in F_2)
                for k in range(j, n):
                    mat[i][k] = (mat[i][k] + mat[j][k]) % 2
                    
    # The matrix is now an identity matrix after elimination, det is 1.
    # If we ever returned 0, it was singular.
    return 1

def solve_loopless_cycle_cover_parity(num_vertices, edges):
    """
    Calculates the parity of the number of loopless cycle covers.
    
    Args:
        num_vertices (int): The number of vertices in the graph (labeled 0 to n-1).
        edges (list of tuples): A list of directed edges, e.g., [(0, 1), (1, 0)].
    """
    if num_vertices == 0:
        print("Graph has no vertices. Number of cycle covers is 1 (the empty cover). Parity is 1.")
        # By convention, det of 0x0 matrix is 1
        return
        
    n = num_vertices
    # Step 1: Create adjacency matrix A over F_2
    A = [[0] * n for _ in range(n)]
    for u, v in edges:
        A[u][v] = 1

    # Step 2: Create matrix M based on the formula
    M = [[0] * n for _ in range(n)]

    # Calculate diagonal elements of M
    for i in range(n):
        diag_sum = 0
        for j in range(n):
            diag_sum = (diag_sum + A[i][j] * A[j][i]) % 2
        M[i][i] = diag_sum

    # Calculate off-diagonal elements of M
    for i in range(n):
        for j in range(n):
            if i != j:
                M[i][j] = (A[i][j] + A[i][j] * A[j][i]) % 2

    print(f"The number of vertices is {n}.")
    print("The adjacency matrix A is:")
    for row in A:
        print(row)
    
    print("\nThe matrix M, whose determinant gives the parity, is:")
    for row in M:
        print(row)
    
    # Step 3: Compute the determinant of M mod 2
    parity = det_mod2(M)
    
    # Step 4: Output the result
    final_eq_str = " * ".join(str(M[i][i]) for i in range(n))
    
    # To demonstrate the determinant calculation, we can show the result on the triangularized form
    # but the provided function directly computes it. For simplicity, we just show the final result.
    # The 'equation' might refer to the product of diagonals of the final triangular matrix.
    # After Gaussian elimination on M, we would get the identity matrix if det=1.
    # The product of its diagonal entries would be 1 * 1 * ... * 1 = 1.
    # If the determinant is 0, the product will be 0.
    if parity == 1:
        print(f"\nThe matrix M is non-singular modulo 2.")
        print(f"The final equation after triangularization would result in det(M) = {' * '.join(['1'] * n)} = 1 (mod 2)")
    else:
        print(f"\nThe matrix M is singular modulo 2.")
        print(f"The final equation after triangularization would have at least one zero on the diagonal, so det(M) = 0 (mod 2)")

    print(f"\nResult: The parity of the number of loopless cycle covers is {parity}.")


if __name__ == '__main__':
    # Example from thought process:
    # V={0,1,2,3}, E={(0,1),(1,0), (2,3),(3,2), (0,2),(1,3)}
    # Vertices are 0-indexed for the code.
    num_v = 4
    edge_list = [(0, 1), (1, 0), (2, 3), (3, 2), (0, 2), (1, 3)]
    
    solve_loopless_cycle_cover_parity(num_v, edge_list)
