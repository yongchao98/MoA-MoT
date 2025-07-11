import numpy as np
import networkx as nx

def analyze_graph(G, name):
    """
    Analyzes a graph to check the condition from the problem and evaluate the answer choices.
    """
    print(f"--- Analyzing graph: {name} ---")
    
    if G.number_of_nodes() <= 3:
        print("Graph has 3 or fewer nodes, but problem requires more.")
        return

    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    # Get connected components and count how many are bipartite
    components = list(nx.connected_components(G))
    c = len(components)
    c_B = 0
    for comp_nodes in components:
        subgraph = G.subgraph(comp_nodes).copy()
        if nx.is_bipartite(subgraph):
            c_B += 1
    
    # Calculate k = dim(null(B^T B)) using the formula k = m - n + c_B
    k = m - n + c_B
    
    # Calculate the largest eigenvalue of the Laplacian
    # L = D - A
    L = nx.laplacian_matrix(G).toarray()
    # Eigenvalues are sorted in ascending order by eigvalsh
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]

    print(f"Graph properties: Nodes n={n}, Edges m={m}, Components c={c}, Bipartite Comp. c_B={c_B}")
    print(f"Calculated k = dim(null(B^T B)) = m - n + c_B = {m} - {n} + {c_B} = {k}")
    print(f"Largest Laplacian eigenvalue λ_n = {lambda_n:.4f}")
    
    # Check if the premise holds: k = λ_n / 2
    premise_holds = np.isclose(k, lambda_n / 2)
    print(f"Premise 'k = λn/2':  {k} = {lambda_n / 2:.4f}  =>  {premise_holds}")
    
    if not premise_holds:
        print("This graph does not satisfy the premise.")
        print("-" * 30 + "\n")
        return
        
    print("\nThis graph SATISFIES the premise. Checking answer choices:")
    # Answer B: The graph has at least k connected components.
    # This translates to c >= k
    print(f"Choice B (c >= k): {c} >= {k} => {c >= k}")
    
    # Answer C: The graph is connected.
    # This translates to c = 1
    print(f"Choice C (c = 1):  {c} = 1 => {c == 1}")
    
    # Answer D: The graph has exactly k connected components.
    # This translates to c = k
    print(f"Choice D (c = k):  {c} = {k} => {c == k}")
    
    print("-" * 30 + "\n")


# 1. A graph that refutes choices B and D
G1 = nx.complete_bipartite_graph(2, 4)
analyze_graph(G1, "K_{2,4}")

# 2. A graph that refutes choice C
G2 = nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))
analyze_graph(G2, "2 disjoint C_4 graphs")
