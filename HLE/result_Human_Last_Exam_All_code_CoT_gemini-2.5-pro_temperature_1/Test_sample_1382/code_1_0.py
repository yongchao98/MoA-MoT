import numpy as np
import networkx as nx

def analyze_graph(G):
    """
    Analyzes a graph to check the relationship given in the problem.
    """
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    if n <= 3:
        print("Graph must have more than 3 nodes.")
        return

    # 1. Calculate null(B^T B) = m - n + k
    # k is the number of connected components
    k = nx.number_connected_components(G)
    nullity_BtB = m - n + k

    # 2. Calculate lambda_n(G)/2
    # Get the graph Laplacian L = D - A
    L = nx.laplacian_matrix(G).toarray()
    # Get the eigenvalues of L
    eigenvalues = np.linalg.eigvalsh(L)
    # Get the largest eigenvalue, lambda_n
    lambda_n = eigenvalues[-1]
    half_lambda_n = lambda_n / 2

    print(f"Graph: {G.graph.get('name', 'Custom Graph')}")
    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (k): {k}")
    print("-" * 20)
    print(f"The equation we want to check is: m - n + k = lambda_n / 2")
    print(f"Left side (m - n + k): {m} - {n} + {k} = {nullity_BtB}")
    print(f"Right side (lambda_n / 2): {lambda_n:.4f} / 2 = {half_lambda_n:.4f}")
    
    if np.isclose(nullity_BtB, half_lambda_n):
        print("\nThe premise holds for this graph.")
    else:
        print("\nThe premise does NOT hold for this graph.")

    # 3. Check the proposed answer B: k >= lambda_n / 2
    print("\nChecking Answer B: k >= lambda_n / 2")
    print(f"Is {k} >= {half_lambda_n:.4f}? {'Yes' if k >= half_lambda_n else 'No'}")

# Create a C4 graph, which we identified as satisfying the premise
C4 = nx.cycle_graph(4)
C4.graph['name'] = 'Cycle graph C4'
analyze_graph(C4)
