import numpy as np
import networkx as nx

def analyze_graph(G, graph_name):
    """
    Analyzes a graph G to check the given condition and the answer choices.
    """
    print(f"--- Analyzing graph: {graph_name} ---")

    # Graph properties
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    # Check if graph has more than 3 nodes
    if n <= 3:
        print(f"Graph has {n} nodes, condition not met.\n")
        return

    is_connected = nx.is_connected(G)
    num_components = nx.number_connected_components(G)
    
    components = list(nx.connected_components(G))
    bipartite_components = 0
    for comp_nodes in components:
        subgraph = G.subgraph(comp_nodes)
        if nx.is_bipartite(subgraph):
            bipartite_components += 1

    # 1. Calculate nullity of B^T * B
    # nullity(B^T * B) = m - n + c_b
    null_B_T_B = m - n + bipartite_components
    print(f"n (nodes) = {n}")
    print(f"m (edges) = {m}")
    print(f"Number of connected components (c) = {num_components}")
    print(f"Number of bipartite components (c_b) = {bipartite_components}")
    print(f"Calculated null(B^T * B) = m - n + c_b = {m} - {n} + {bipartite_components} = {null_B_T_B}")
    
    # 2. Calculate largest eigenvalue of the Laplacian
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]
    
    print(f"Largest Laplacian eigenvalue (lambda_n) = {lambda_n:.4f}")
    
    # 3. Check if the condition from the problem holds
    # null(B^T * B) == lambda_n / 2
    condition_holds = np.isclose(null_B_T_B, lambda_n / 2)
    print(f"Condition check: {null_B_T_B} == {lambda_n / 2:.4f}? {'Yes' if condition_holds else 'No'}")
    
    if not condition_holds:
        print("This graph does not satisfy the condition.\n")
        return
        
    print("\nThis graph SATISFIES the condition. Checking answer choices:")
    # Check Choice B: The graph has at least lambda_n/2 connected components
    choice_b = num_components >= lambda_n / 2
    print(f"Choice B (c >= lambda_n/2): {num_components} >= {lambda_n/2:.2f}? {choice_b}")

    # Check Choice C: The graph is connected
    choice_c = is_connected
    print(f"Choice C (is connected): {choice_c}")
    
    # Check Choice D: The graph has exactly lambda_n/2 connected components
    choice_d = np.isclose(num_components, lambda_n / 2)
    print(f"Choice D (c == lambda_n/2): {num_components} == {lambda_n/2:.2f}? {choice_d}\n")


# Case 1: K4 (Complete graph on 4 vertices)
G1 = nx.complete_graph(4)
analyze_graph(G1, "K4")

# Case 2: K4 U K1 (Disconnected graph)
G2 = nx.disjoint_union(nx.complete_graph(4), nx.empty_graph(1))
analyze_graph(G2, "K4 U K1")
