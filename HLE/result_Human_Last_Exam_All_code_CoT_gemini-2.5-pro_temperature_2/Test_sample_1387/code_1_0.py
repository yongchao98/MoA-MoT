import numpy as np

def solve_loopless_cycle_cover_parity(adj_matrix):
    """
    Calculates the parity of the number of loopless cycle covers of a directed graph.

    Args:
        adj_matrix (list of lists or numpy.ndarray): The adjacency matrix of the graph.
    """
    # Ensure we are using a numpy array for easier operations
    A = np.array(adj_matrix, dtype=int)
    n = A.shape[0]

    if A.shape != (n, n):
        raise ValueError("Adjacency matrix must be square.")

    # All computations are modulo 2
    A = A % 2

    # Step 1: Construct the matrix S representing symmetric edges
    # S[i][j] = 1 iff A[i][j] = 1 and A[j][i] = 1
    # This is the Hadamard product of A and its transpose.
    S = (A * A.T) % 2
    
    # Step 2: Construct the matrix B = A + S (mod 2)
    B = (A + S) % 2
    
    print("Given the adjacency matrix A:")
    print(A)
    print("\nMatrix of symmetric edges S = A * A^T (element-wise):")
    print(S)
    print("\nMatrix for determinant calculation B = A + S (mod 2):")
    print(B)

    # Step 3: Compute the determinant of B modulo 2
    # For a matrix over F_2, the determinant can be found via Gaussian elimination.
    # The determinant will be 1 if the matrix is invertible, 0 otherwise.
    
    M = B.copy() # Make a copy for Gaussian elimination
    
    det = 1
    for i in range(n):
        # Find pivot
        pivot = i
        while pivot < n and M[pivot][i] == 0:
            pivot += 1
        
        if pivot == n: # No pivot found in this column
            det = 0
            break
        
        # Swap rows to bring pivot to diagonal
        if pivot != i:
            M[[i, pivot]] = M[[pivot, i]]
            # Swapping rows flips the sign of the determinant, but mod 2, -1=1 so it doesn't matter.

        # Eliminate other 1s in the same column
        for j in range(i + 1, n):
            if M[j][i] == 1:
                M[j, :] = (M[j, :] + M[i, :]) % 2

    print(f"\nThe determinant of B is {det}.")
    print("\nThe parity of the number of loopless cycle covers is (0 for even, 1 for odd).")
    print("Final result:")
    print(det)

if __name__ == '__main__':
    # Example usage with a simple directed graph
    # Let's consider a graph with 4 vertices.
    # Edges: 0->1, 1->0 (a 2-cycle), 1->2, 2->3, 3->0, 3->1
    example_adj_matrix = [
        [0, 1, 0, 0],
        [1, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 1, 0, 0]
    ]
    solve_loopless_cycle_cover_parity(example_adj_matrix)
