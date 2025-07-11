import numpy as np

def solve():
    """
    This function verifies the premise for the graph K4 and prints the calculations.
    """
    # Define the graph K4
    n = 4  # number of nodes
    m = 6  # number of edges
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    # 1. Construct the unoriented incidence matrix B
    B = np.zeros((n, m))
    for i, edge in enumerate(edges):
        B[edge[0], i] = 1
        B[edge[1], i] = 1

    # 2. Calculate nullity of B
    # Nullity = number of columns - rank
    rank_B = np.linalg.matrix_rank(B)
    nullity_B = m - rank_B

    # 3. Construct the Laplacian L
    A = np.zeros((n, n))
    for i, j in edges:
        A[i, j] = A[j, i] = 1
    D = np.diag(np.sum(A, axis=1))
    L = D - A

    # 4. Find the largest eigenvalue of L
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]

    # 5. Check the premise and print the results
    lhs = nullity_B
    rhs = lambda_n / 2

    print("Someone hands you a graph and tells you that the following equation holds:")
    print("null(B^T B) = lambda_n(G) / 2")
    print("\nLet's test this premise for the complete graph K4.")
    print(f"The number of nodes n = {n}")
    print(f"The number of edges m = {m}")
    print(f"The number of bipartite components k_bip = 0")
    print(f"Based on the formula nullity(B) = m - n + k_bip, the expected nullity is {m - n + 0}.")
    print("\nVerifying with code:")
    print(f"Calculated nullity of B (LHS of the equation): {lhs}")
    print(f"Largest eigenvalue of the Laplacian, lambda_n: {lambda_n:.4f}")
    print(f"Calculated lambda_n / 2 (RHS of the equation): {rhs}")

    if np.isclose(lhs, rhs):
        print("\nThe premise holds for K4.")
    else:
        print("\nThe premise does NOT hold for K4.")
        
    print("\nFor K4, which is connected, the answer choice 'C. The graph is connected' is TRUE.")
    print("The other options are false for K4.")

solve()