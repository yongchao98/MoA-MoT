import networkx as nx
import numpy as np

def analyze_graph_property(G, graph_name):
    """
    Analyzes a graph to test the property given in the problem.
    The property is null(B^T B) = lambda_n(G) / 2, which translates to
    m - n + k_b = lambda_n(G) / 2.

    We test this against Answer Choice B: k >= lambda_n(G) / 2
    """
    print(f"--- Analyzing graph: {graph_name} ---")

    # 1. Get graph invariants
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    # Get connected components and count them
    components = list(nx.connected_components(G))
    k = len(components)
    
    # Determine the number of bipartite components
    k_b = 0
    for i, comp_nodes in enumerate(components):
        subgraph = G.subgraph(comp_nodes)
        if nx.is_bipartite(subgraph):
            k_b += 1
    
    # 2. Calculate nullity(B^T B) using the formula m - n + k_b
    # This is the left-hand side (LHS) of the condition.
    nullity_B_T_B = m - n + k_b

    # 3. Calculate the largest eigenvalue of the Laplacian, lambda_n(G)
    # This is needed for the right-hand side (RHS) of the condition.
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    # Sort eigenvalues to find the largest, lambda_n
    eigenvalues.sort()
    lambda_n = eigenvalues[-1]
    
    rhs = lambda_n / 2

    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (k): {k}")
    print(f"Number of bipartite components (k_b): {k_b}")
    print(f"Largest Laplacian eigenvalue (lambda_n): {lambda_n:.4f}")
    
    print("\nChecking if the graph satisfies the given condition...")
    print(f"Condition is: null(B^T B) = lambda_n / 2")
    print(f"Which translates to: m - n + k_b = lambda_n / 2")
    print(f"LHS (m - n + k_b) = {m} - {n} + {k_b} = {nullity_B_T_B}")
    print(f"RHS (lambda_n / 2) = {lambda_n:.4f} / 2 = {rhs:.4f}")
    
    # Use a tolerance for float comparison
    if np.isclose(nullity_B_T_B, rhs):
        print("Conclusion: The graph satisfies the condition.")
        
        # 4. Now test Answer Choice B
        print("\nTesting Answer Choice B: The graph has at least lambda_n/2 connected components.")
        print(f"This means we must check if k >= lambda_n / 2")
        print(f"Check: {k} >= {rhs:.4f}?")
        if k >= rhs:
            print("Result: The inequality holds for this graph.")
        else:
            print("Result: The inequality DOES NOT hold for this graph.")
            print(f"Therefore, {graph_name} is a counterexample to Answer Choice B.")
    else:
        print("Conclusion: The graph DOES NOT satisfy the condition, so it cannot be used as a counterexample.")

# Let's use the complete graph on 4 nodes, K_4, as a potential counterexample.
# K_4 has n=4 > 3, so it's a valid graph.
G_K4 = nx.complete_graph(4)
analyze_graph_property(G_K4, "K_4 (Complete Graph on 4 nodes)")

print("\n--- Final Analysis ---")
print("The analysis for K_4 shows that a graph can satisfy the given condition null(B^T B) = lambda_n/2.")
print("However, for K_4, the number of connected components (k=1) is NOT greater than or equal to lambda_n/2 (which is 2).")
print("This disproves option B.")
print("Similar counterexamples can be constructed to disprove options A, C, and D.")
print("Therefore, none of the options A, B, C, or D are correct.")

<<<E>>>