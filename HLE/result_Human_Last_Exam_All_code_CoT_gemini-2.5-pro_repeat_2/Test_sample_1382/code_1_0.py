import numpy as np

def get_k4_properties():
    """
    Calculates properties for the complete graph K4.
    """
    # Properties for K4
    n = 4  # number of nodes
    m = 6  # number of edges
    c = 1  # number of connected components

    # Adjacency matrix for K4
    A = np.ones((n, n)) - np.identity(n)

    # Degree matrix for K4
    D = np.diag(np.sum(A, axis=1))

    # Laplacian matrix L = D - A
    L = D - A

    # Eigenvalues of the Laplacian
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = max(eigenvalues)

    # Cyclomatic number (beta)
    beta = m - n + c
    
    return beta, lambda_n

def check_inequality():
    """
    Checks if beta > lambda_n / 2 for K4.
    """
    beta, lambda_n = get_k4_properties()
    
    lhs = beta
    rhs = lambda_n / 2
    
    print(f"For the complete graph K4:")
    print(f"Number of nodes n = 4")
    print(f"Number of edges m = 6")
    print(f"Number of connected components c = 1")
    print(f"Cyclomatic number beta = m - n + c = {6} - {4} + {1} = {beta}")
    print(f"Largest Laplacian eigenvalue lambda_n = {lambda_n:.4f}")
    print("\nChecking the inequality 2 * beta > lambda_n, which is beta > lambda_n / 2:")
    print(f"Left side (beta): {lhs}")
    print(f"Right side (lambda_n / 2): {rhs}")
    
    if lhs > rhs:
        print(f"Result: {lhs} > {rhs} is TRUE.")
        print("\nThis demonstrates a case where a component cannot satisfy the condition required for a disconnected graph,")
        print("supporting the conclusion that the graph must be connected.")
    else:
        print(f"Result: {lhs} > {rhs} is FALSE.")

check_inequality()
