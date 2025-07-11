import numpy as np

def inv_gf2(A):
    """
    Computes the inverse of a matrix over GF(2).
    Returns the inverse matrix or None if the matrix is singular.
    """
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise ValueError("Input matrix must be square")
    
    A_copy = A.astype(int)
    I = np.eye(n, dtype=int)
    A_ext = np.hstack([A_copy, I])

    for i in range(n):
        # Find a pivot in column i at or below row i
        pivot_row = i
        while pivot_row < n and A_ext[pivot_row, i] == 0:
            pivot_row += 1

        if pivot_row == n:
            # No pivot found, matrix is singular
            return None

        # Swap the current row with the pivot row
        A_ext[[i, pivot_row]] = A_ext[[pivot_row, i]]

        # Eliminate other 1s in the current column i
        for j in range(n):
            if i != j and A_ext[j, i] == 1:
                A_ext[j, :] = (A_ext[j, :] + A_ext[i, :]) % 2

    # The right part of A_ext is now the inverse
    return A_ext[:, n:]

def det_gf2(A):
    """
    Computes the determinant of a matrix over GF(2).
    """
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise ValueError("Input matrix must be square")
        
    M = A.copy().astype(int)

    for i in range(n):
        # Find a pivot in column i at or below row i
        pivot_row = i
        while pivot_row < n and M[pivot_row, i] == 0:
            pivot_row += 1

        if pivot_row == n:
            # No pivot in this column, so the determinant is 0
            return 0

        # Swap rows to bring pivot to the diagonal
        M[[i, pivot_row]] = M[[pivot_row, i]]

        # Use the pivot to eliminate other 1s below it in the same column
        for j in range(i + 1, n):
            if M[j, i] == 1:
                M[j, :] = (M[j, :] + M[i, :]) % 2

    # After Gaussian elimination, the matrix is upper triangular.
    # The determinant is the product of the diagonal elements.
    # In our case, all diagonal elements are 1, so the determinant is 1.
    return 1

def solve_loopless_cycle_cover_parity(G, n):
    """
    Solves the LooplessCycleCover parity problem.
    G: list of arcs (u, v) representing the directed multigraph.
    n: number of vertices (labeled 0 to n-1).
    """
    print(f"Solving for a graph with {n} vertices.")
    
    # 1. Create adjacency matrix M and reduce it mod 2 to get A
    M = np.zeros((n, n), dtype=int)
    for u, v in G:
        M[u, v] += 1
    A = M % 2
    print("\nAdjacency matrix mod 2 (A):\n", A)

    # 2. Compute S where S_ij = A_ij * A_ji
    S = A * A.T
    print("\nMatrix of 2-cycle pairs (S):\n", S)
    
    # 3. Compute B = I + S
    I = np.eye(n, dtype=int)
    B = (I + S) % 2
    print("\nMatrix B = I + S (mod 2):\n", B)

    # 4. Compute inverse of B
    B_inv = inv_gf2(B)
    
    if B_inv is None:
        print("\nMatrix B is singular. The number of loopless cycle covers is even.")
        print("Final Result: The parity of loopless cycle covers is 0.")
        return 0

    print("\nInverse of B (mod 2):\n", B_inv)

    # 5. Compute T = B_inv @ A
    T = (B_inv @ A) % 2
    print("\nMatrix T = (B_inv @ A) (mod 2):\n", T)
    
    # 6. Compute determinant of T
    result = det_gf2(T)
    
    print(f"\nThe determinant of T is {result}.")
    print(f"Final equation: parity = det(T) mod 2 = {result}")
    print(f"\nFinal Result: The parity of loopless cycle covers is: {result}")
    return result

if __name__ == '__main__':
    # Define an example directed graph
    # V = {0, 1, 2, 3}, n = 4
    # E is a list of arcs, can be a multiset
    example_graph_edges = [
        (0, 1), (1, 0),  # A 2-cycle
        (1, 2), (2, 1),  # Another 2-cycle
        (2, 3),
        (3, 0),
        (0, 2)
    ]
    num_vertices = 4
    
    # Let's manually check the loopless cycle covers for this example.
    # A loopless cycle cover must be composed of cycles of length >= 3.
    # For n=4, the only possibility is a single 4-cycle.
    # Potential 4-cycles (permutations of (0,1,2,3)):
    # 1. 0 -> 1 -> 2 -> 3 -> 0. Arcs: (0,1), (1,2), (2,3), (3,0). All exist. This is one valid cover.
    # Any other 4-cycle, like 0->2->3->1->0 involves arc (3,1), which doesn't exist.
    # There is exactly 1 loopless cycle cover. The parity is 1.
    
    solve_loopless_cycle_cover_parity(example_graph_edges, num_vertices)
