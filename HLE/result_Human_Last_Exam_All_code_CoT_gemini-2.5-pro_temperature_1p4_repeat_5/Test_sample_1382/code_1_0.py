import networkx as nx
import numpy as np

def analyze_graph(G, name):
    """
    Analyzes a given graph according to the problem statement and checks the answer choices.
    """
    print(f"--- Analyzing graph: {name} ---")

    n = G.number_of_nodes()
    m = G.number_of_edges()

    if n <= 3:
        print("Graph has 3 or fewer nodes, skipping as per problem statement.\n")
        return

    # Get connected components
    components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    k = len(components)

    # Calculate k_bipartite
    k_bipartite = 0
    for comp in components:
        if nx.is_bipartite(comp):
            k_bipartite += 1
    
    # Calculate lambda_n (largest eigenvalue of Laplacian)
    # Note: networkx returns them sorted ascendingly
    try:
        laplacian_eigenvalues = nx.laplacian_spectrum(G)
        lambda_n = laplacian_eigenvalues[-1]
    except nx.NetworkXError:
        print("Could not compute Laplacian spectrum (graph might be empty).")
        return


    # Check the given property: null(B^T B) = lambda_n / 2
    # which is equivalent to m - n + k_bipartite = lambda_n / 2
    lhs = m - n + k_bipartite
    rhs = lambda_n / 2.0
    
    print(f"Parameters: n={n}, m={m}, k={k}, k_bipartite={k_bipartite}")
    print(f"Largest Laplacian eigenvalue (lambda_n) = {lambda_n:.4f}")
    print(f"Checking property: m - n + k_bipartite = {lhs}")
    print(f"lambda_n / 2 = {rhs:.4f}")
    
    # Use a tolerance for float comparison
    property_holds = np.isclose(lhs, rhs)
    print(f"Does the property hold for this graph? {'Yes' if property_holds else 'No'}")

    if not property_holds:
        print("Property does not hold, so this graph is not a valid test case.\n")
        return

    print("\nThis graph satisfies the property. Now checking the answer choices:")
    
    # Check C: The graph is connected
    c_holds = (k == 1)
    print(f"Choice C (is connected?): {c_holds}")
    
    # Check B: The graph has at least lambda_n/2 connected components
    b_holds = (k >= rhs)
    print(f"Choice B (k >= lambda_n/2?): {k} >= {rhs:.4f} -> {b_holds}")

    # Check D: The graph has exactly lambda_n/2 connected components
    d_holds = np.isclose(k, rhs)
    print(f"Choice D (k == lambda_n/2?): {k} == {rhs:.4f} -> {d_holds}")
    
    print("-" * 25 + "\n")

# --- Main Execution ---

# Graph 1: K_4 (Complete graph on 4 nodes)
# It is connected and not bipartite.
G1 = nx.complete_graph(4)
analyze_graph(G1, "K_4")

# Graph 2: K_4 U P_2 (A disconnected graph)
# One component is K_4, the other is a single edge (P_2).
G2 = nx.disjoint_union(nx.complete_graph(4), nx.path_graph(2))
analyze_graph(G2, "K_4 U P_2")

# Manual check for Choice A on K_4:
# For K_4, lambda_n/2 = 2. We need to drop 2 edges.
# K_4 has 6 edges. If we drop 2 non-adjacent edges, we are left with a cycle C_4.
# In a C_4, every node has degree 2. So, no node has degree <= 1.
# Therefore, Choice A is false for K_4.

print("--- Summary of Findings ---")
print("1. K_4 satisfies the property. For this graph, choices B and D are FALSE. As shown manually, choice A is also FALSE.")
print("2. K_4 U P_2 satisfies the property. For this graph, choice C is FALSE.")
print("\nSince no single choice (A, B, C, or D) is true for all graphs that satisfy the given condition, none of them is a necessary consequence.")