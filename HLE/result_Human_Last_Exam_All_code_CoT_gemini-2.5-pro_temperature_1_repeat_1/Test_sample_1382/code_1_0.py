import networkx as nx
import numpy as np

def analyze_graph(G, graph_name):
    """
    Analyzes a graph to check if it satisfies the given condition and the answer choices.
    """
    print(f"--- Analyzing graph: {graph_name} ---")

    # Graph properties
    n = G.number_of_nodes()
    m = G.number_of_edges()
    k = nx.number_connected_components(G)
    
    # Calculate k_bipartite
    k_bipartite = 0
    if k > 0:
        for component_nodes in nx.connected_components(G):
            subgraph = G.subgraph(component_nodes)
            if nx.is_bipartite(subgraph):
                k_bipartite += 1
    
    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (k): {k}")
    print(f"Number of bipartite components (k_bipartite): {k_bipartite}")

    # Calculate nullity of B^T * B
    # null(B^T * B) = m - n + k_bipartite
    nullity_btb = m - n + k_bipartite
    print(f"null(B^T * B) = m - n + k_bipartite = {m} - {n} + {k_bipartite} = {nullity_btb}")
    
    # Calculate largest eigenvalue of the Laplacian
    # Note: np.linalg.eigvalsh is efficient for symmetric matrices
    # networkx returns a sorted list of eigenvalues
    laplacian_eigenvalues = nx.laplacian_spectrum(G)
    lambda_n = laplacian_eigenvalues[-1]
    # Use a small tolerance for float comparison
    lambda_n = round(lambda_n, 8)
    print(f"Largest Laplacian eigenvalue (lambda_n): {lambda_n}")
    
    # Check the condition from the problem
    val_from_condition = lambda_n / 2
    print(f"Value from condition (lambda_n / 2): {val_from_condition}")
    
    if np.isclose(nullity_btb, val_from_condition):
        print("\nThe graph satisfies the condition null(B^T * B) = lambda_n / 2.")
        
        # Now check the answer choices
        print("\n--- Checking Answer Choices ---")
        
        # Choice B: The graph has at least lambda_n/2 connected components
        choice_b_holds = k >= val_from_condition
        print(f"Choice B (k >= lambda_n/2): Is {k} >= {val_from_condition}? {choice_b_holds}")
        
        # Choice D: The graph has exactly lambda_n/2 connected components
        choice_d_holds = np.isclose(k, val_from_condition)
        print(f"Choice D (k == lambda_n/2): Is {k} == {val_from_condition}? {choice_d_holds}")

    else:
        print("The graph does not satisfy the condition.")

# Create the complete graph K4
G_k4 = nx.complete_graph(4)

# Run the analysis
analyze_graph(G_k4, "K4 (Complete Graph on 4 nodes)")
