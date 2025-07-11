import numpy as np

def analyze_graph():
    """
    Analyzes the graph G = C_4 U C_4 to verify the problem statement.
    """
    # Define graph properties
    n = 8  # Number of nodes
    m = 8  # Number of edges
    k = 2  # Number of connected components

    # Adjacency matrix for two disjoint C4 cycles
    # Nodes 0-3 are the first C4, nodes 4-7 are the second C4
    A = np.array([
        # First C4
        [0, 1, 0, 1, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0],
        # Second C4
        [0, 0, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 1, 0, 1, 0]
    ])

    # Calculate Degree and Laplacian matrices
    D = np.diag(np.sum(A, axis=1))
    L = D - A

    # Calculate eigenvalues of the Laplacian
    # Use eigvalsh for symmetric matrices for better performance and precision
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = np.max(eigenvalues)

    # The nullity of B^T B is the cyclomatic number mu
    mu = m - n + k
    
    # The value from the right side of the equation
    rhs_value = lambda_n / 2

    print(f"Graph properties for G = C4 U C4:")
    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (k): {k}")
    print("-" * 30)

    print(f"Equation from the problem: null(B^T B) = lambda_n(G) / 2")
    print(f"LHS (null(B^T B) = m - n + k): {m} - {n} + {k} = {mu}")
    print(f"RHS (lambda_n(G) / 2): {lambda_n:.2f} / 2 = {rhs_value:.2f}")

    if np.isclose(mu, rhs_value):
        print("\nThe given equation holds true for this graph.")
    else:
        print("\nThe given equation DOES NOT hold true for this graph.")

    print("\nChecking Answer Choice D: The graph has exactly lambda_n(G) / 2 connected components.")
    print(f"Number of components k = {k}")
    print(f"Value of lambda_n(G) / 2 = {rhs_value:.2f}")
    if np.isclose(k, rhs_value):
        print("Answer choice D is consistent with this example.")
    else:
        print("Answer choice D is NOT consistent with this example.")

analyze_graph()