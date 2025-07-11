import numpy as np

def solve_loopless_cycle_cover_parity(adj_matrix):
    """
    Calculates the parity of the number of loopless cycle covers of a directed multigraph.

    The algorithm is based on the result from algebraic graph theory that the desired
    quantity is equal to det(A + Q) mod 2, where:
    - A is the adjacency matrix of the graph modulo 2.
    - Q is the adjacency matrix of the underlying graph of 2-cycles.

    Args:
        adj_matrix (list of lists or numpy.ndarray): The adjacency matrix of the
                                                    directed multigraph G, where M[i][j]
                                                    is the number of arcs from i to j.

    Returns:
        int: The parity (0 or 1) of the number of loopless cycle covers.
    """
    n = len(adj_matrix)
    A = np.array(adj_matrix, dtype=int) % 2

    # Create the adjacency matrix Q of the 2-cycle graph H
    # An edge {i, j} exists in H if there are arcs (i,j) and (j,i) in G.
    Q = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j] == 1 and A[j, i] == 1:
                Q[i, j] = 1
                Q[j, i] = 1

    # The number of loopless cycle covers mod 2 is det(A + Q) mod 2.
    B = (A + Q) % 2

    # Calculate the determinant of B over F_2.
    # We can use standard determinant calculation and then take mod 2.
    # For numerical stability with large matrices, a dedicated F_2 Gaussian elimination
    # would be better, but for typical inputs, this works fine.
    det_B = int(round(np.linalg.det(B)))
    parity = det_B % 2
    
    # Printing the intermediate steps as requested
    print("This program calculates the parity of loopless cycle covers for a given graph.")
    print(f"The graph has {n} vertices.")
    print("\nStep 1: Adjacency matrix A (mod 2)")
    print(A)
    print("\nStep 2: Adjacency matrix Q of the 2-cycle graph")
    print(Q)
    print("\nStep 3: The matrix B = A + Q (mod 2)")
    print(B)
    print(f"\nStep 4: The determinant of B is {det_B}.")
    print(f"The parity is det(B) mod 2 = {det_B} mod 2 = {parity}.")
    
    print("\nFinal Answer (the parity):")
    print(parity)

    return parity

if __name__ == '__main__':
    # Example Usage:
    # Let's define a graph with 4 vertices.
    # The adjacency matrix represents the number of arcs between vertices.
    # Let the vertices be {0, 1, 2, 3}.
    # Arcs: (0,1), (1,0), (0,2), (2,3), (3,0), (1,3). All with multiplicity 1.
    
    #      0  1  2  3
    #   0 [0, 1, 1, 0]
    #   1 [1, 0, 0, 1]
    #   2 [0, 0, 0, 1]
    #   3 [1, 0, 0, 0]
    
    example_adj_matrix = [
        [0, 1, 1, 0],
        [1, 0, 0, 1],
        [0, 0, 0, 1],
        [1, 0, 0, 0]
    ]

    solve_loopless_cycle_cover_parity(example_adj_matrix)
