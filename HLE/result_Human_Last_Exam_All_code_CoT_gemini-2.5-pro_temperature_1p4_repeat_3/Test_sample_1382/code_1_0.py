import numpy as np

def analyze_graph_properties(edges, num_nodes):
    """
    Analyzes a graph to compute properties from the problem description.

    Args:
        edges (list of tuples): A list of edges, e.g., [(0, 1), (1, 2)].
        num_nodes (int): The total number of nodes in the graph.
    """
    if num_nodes == 0:
        print("Graph has no nodes.")
        return

    num_edges = len(edges)
    
    # 1. Construct the unoriented incidence matrix B (n x m)
    B = np.zeros((num_nodes, num_edges), dtype=int)
    for i, edge in enumerate(edges):
        u, v = edge
        if u < num_nodes and v < num_nodes:
            B[u, i] = 1
            B[v, i] = 1
    
    # 2. Construct Adjacency and Degree matrices to get the Laplacian L
    A = np.zeros((num_nodes, num_nodes), dtype=int)
    for u, v in edges:
         if u < num_nodes and v < num_nodes:
            A[u, v] = 1
            A[v, u] = 1
            
    D = np.diag(np.sum(A, axis=1))
    L = D - A

    # 3. Calculate null(L) = k (number of connected components)
    # The nullity is the number of eigenvalues that are close to zero.
    eigenvalues_L = np.linalg.eigvalsh(L)
    k = np.sum(np.isclose(eigenvalues_L, 0))

    # 4. Calculate lambda_n(G) (largest eigenvalue of L)
    lambda_n = eigenvalues_L[-1] if len(eigenvalues_L) > 0 else 0

    # 5. Calculate null(B^T B)
    # We can do this without computing m-n+b by finding the nullity directly.
    if num_edges > 0:
        B_transpose_B = B.T @ B
        eigenvalues_BTB = np.linalg.eigvalsh(B_transpose_B)
        null_B_T_B = np.sum(np.isclose(eigenvalues_BTB, 0))
    else:
        null_B_T_B = 0

    print("--- Graph Analysis ---")
    print(f"Number of nodes (n): {num_nodes}")
    print(f"Number of edges (m): {num_edges}")
    print("\nInterpreting the statement:")

    # Based on our reasoning, the statement implies k = lambda_n / 2
    # The code below prints out the values for you to see.
    # Note: A random graph will not satisfy the condition, the problem says
    # "Someone hands you a graph" where this condition is true.
    
    print("\n1. Number of connected components (k):")
    print(f"The number of connected components is k = null(L) = {k}")

    print("\n2. The value from the right side of the equation:")
    print(f"The largest eigenvalue of the Laplacian is lambda_n = {lambda_n:.4f}")
    if lambda_n > 0:
        print(f"So, the value lambda_n / 2 = {lambda_n/2:.4f}")
    
    print("\nIf the statement k = lambda_n / 2 were true for this graph:")
    # We can now print the final equation k = lambda_n / 2 with the numbers for our example graph
    print(f"The equation would be: {k} = {lambda_n:.4f} / 2")
    
    print("\n--- For context, the literal value of null(B^T B) ---")
    print(f"The actual nullity of B^T B for this graph is {null_B_T_B}")
    

# Example Graph: Two disconnected edges (P2 U P2).
# Nodes 0-1 are an edge, and nodes 2-3 are an edge.
# n = 4, k = 2
example_edges = [(0, 1), (2, 3)]
example_nodes = 4
analyze_graph_properties(example_edges, example_nodes)

<<<D>>>