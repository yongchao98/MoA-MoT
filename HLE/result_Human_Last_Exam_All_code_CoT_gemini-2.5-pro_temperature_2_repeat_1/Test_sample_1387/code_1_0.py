import numpy as np

def solve_loopless_cycle_cover_parity():
    """
    Calculates the parity of the number of loopless cycle covers for a given directed graph.

    The problem is solved by constructing a specific matrix and calculating its determinant
    over GF(2). The parity of loopless cycle covers in a directed graph G is equal to
    the parity of the number of matchings in its underlying graph of symmetric edges S_G.
    This, in turn, is given by det(I + A_S) mod 2, where A_S is the adjacency matrix of S_G.
    """

    # Step 1: Define the directed multigraph G.
    # We use an adjacency matrix where M[i, j] is the number of arcs from i+1 to j+1.
    # Example Graph G has vertices {1, 2, 3, 4} and arcs:
    # (1,2), (2,1), (1,3), (3,4), (4,1), (2,4), (4,2)
    # n is the number of vertices.
    n = 4
    # The adjacency matrix for the number of edges.
    M = np.array([
        [0, 1, 1, 0], # Arcs from vertex 1
        [1, 0, 0, 1], # Arcs from vertex 2
        [0, 0, 0, 1], # Arcs from vertex 3
        [1, 1, 0, 0]  # Arcs from vertex 4
    ])
    
    print("This program computes the parity of the number of loopless cycle covers.")
    print(f"The input graph G has {n} vertices.")
    print("The number of arcs between vertices is given by the matrix M:")
    print(M)
    print("-" * 20)

    # Step 2: Get the adjacency matrix A over GF(2).
    A = M % 2
    print("Step 1: Construct the adjacency matrix A over GF(2) (A[i,j] = #arcs(i,j) mod 2).")
    print("Matrix A:")
    print(A)
    print("-" * 20)

    # Step 3: Construct the adjacency matrix A_S of the graph of symmetric edges S_G.
    # An edge {i,j} exists in S_G if arcs (i,j) and (j,i) both exist in G.
    # A_S[i,j] = 1 iff A[i,j]=1 and A[j,i]=1.
    A_S = A * A.T
    print("Step 2: Construct the adjacency matrix A_S for the graph of symmetric edges S_G.")
    print("A_S = A AND A^T.")
    print("Matrix A_S:")
    print(A_S)
    print("-" * 20)

    # Step 4: Construct the matrix B = I + A_S.
    I = np.identity(n, dtype=int)
    B = (I + A_S) % 2
    print("Step 3: Construct the matrix B = I + A_S, whose determinant mod 2 gives the answer.")
    print("Matrix B:")
    print(B)
    print("-" * 20)
    
    # Step 5: Compute the determinant of B over GF(2).
    # For integer matrices, det(B) is an integer. We take the absolute value and then mod 2.
    det_B = int(round(np.linalg.det(B)))
    parity = abs(det_B) % 2
    
    print("Step 4: Compute the determinant of B. The result modulo 2 is the parity.")
    print(f"The determinant of B is {det_B}.")
    # Here, we show each number in the final equation as requested
    print(f"The final equation is: parity = abs({det_B}) mod 2")
    print(f"{abs(det_B)} mod 2 = {parity}")
    
    print("\nResult:")
    if parity == 1:
        print("The number of loopless cycle covers is ODD.")
    else:
        print("The number of loopless cycle covers is EVEN.")

if __name__ == "__main__":
    solve_loopless_cycle_cover_parity()
