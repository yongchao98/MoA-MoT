import numpy as np

def solve_loopless_cycle_cover_parity(adj_matrix):
    """
    Computes the parity of the number of loopless cycle covers of a directed multigraph G.

    The algorithm is based on an inclusion-exclusion argument and a Pfaffian formula
    for sums of determinants over matchings.

    Args:
        adj_matrix (list of lists or numpy.ndarray): The adjacency matrix of the
                    multigraph G, where adj_matrix[i][j] is the number of arcs
                    from vertex i to vertex j.

    Returns:
        int: The parity (0 or 1) of the number of loopless cycle covers.
    """
    A = np.array(adj_matrix, dtype=int)
    n = A.shape[0]

    # 1. Get the adjacency matrix B over GF(2)
    B = A % 2

    # 2. Construct the adjacency matrix of the undirected graph G''
    # An edge {u,v} exists in G'' if arcs (u,v) and (v,u) both exist in G (in odd numbers).
    A_G_double_prime = (B * B.T)
    
    # 3. Construct the matrix X = B + B^T (mod 2)
    # This is the skew-symmetric part of B.
    X = (B + B.T) % 2
    
    # 4. Construct the top-left block of the final matrix K
    # J = A_G'' ⊙ (B + B^T)
    J = A_G_double_prime * X
    
    # 5. Construct the 2n x 2n matrix K
    K = np.zeros((2 * n, 2 * n), dtype=int)
    I = np.identity(n, dtype=int)
    
    # K = [[J, I],
    #      [I, 0]]
    K[0:n, 0:n] = J
    K[0:n, n:2*n] = I
    K[n:2*n, 0:n] = I
    
    # 6. Compute the determinant of K over GF(2).
    # We use Gaussian elimination for this.
    det_K = 1
    matrix = K.copy()
    
    for i in range(2 * n):
        # Find pivot
        pivot_row = i
        while pivot_row < 2 * n and matrix[pivot_row, i] == 0:
            pivot_row += 1
            
        if pivot_row == 2 * n:
            # No pivot in this column
            det_K = 0
            break
        
        # Swap rows to bring pivot to diagonal
        if pivot_row != i:
            matrix[[i, pivot_row]] = matrix[[pivot_row, i]]
            # Swapping rows changes the sign of determinant, but in GF(2) this doesn't matter.
            
        # det_K is multiplied by the pivot, which is 1.
        
        # Eliminate other 1s in the same column
        for j in range(2 * n):
            if i != j and matrix[j, i] == 1:
                matrix[j, :] = (matrix[j, :] + matrix[i, :]) % 2

    # The determinant is the product of diagonal elements, which are all 1s if non-singular.
    
    final_parity = det_K % 2
    
    print(f"The input graph has {n} vertices.")
    print(f"The parity of the number of loopless cycle covers is: {final_parity}")
    return final_parity

if __name__ == '__main__':
    # Example Usage:
    # A complete graph K4, where every pair of vertices has an arc in both directions.
    # A = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]].
    # Cycle covers: 
    # - 3!/3=2 4-cycles: (1234), (1432). These are loopless.
    # - 3 permutations of type (12)(34). These are not loopless. (12)(34), (13)(24), (14)(23).
    # Total loopless covers = 2. Parity = 0.
    # Let's test this with the code.
    k4_adj = [[0, 1, 1, 1],
              [1, 0, 1, 1],
              [1, 1, 0, 1],
              [1, 1, 1, 0]]
    solve_loopless_cycle_cover_parity(k4_adj)

    print("-" * 20)

    # Example 2: Directed 3-cycle 1->2->3->1
    # A = [[0,1,0],[0,0,1],[1,0,0]]
    # Only one cycle cover: (123). It's loopless. Parity = 1.
    c3_adj = [[0, 1, 0],
              [0, 0, 1],
              [1, 0, 0]]
    solve_loopless_cycle_cover_parity(c3_adj)
    
    print("-" * 20)
    
    # Example 3: A graph with no loopless cycle covers.
    # Two disjoint arcs: 1->2, 3->4. No cycle cover exists. Parity should be 0.
    # A = [[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]]
    no_cc_adj = [[0, 1, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 0, 0]]
    solve_loopless_cycle_cover_parity(no_cc_adj)
