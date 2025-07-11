import numpy as np
import networkx as nx

def analyze_graph_property(G, graph_name):
    """
    Analyzes if a given graph G satisfies the condition and checks the answer choices.
    """
    print(f"--- Analyzing Graph: {graph_name} ---")
    
    n = G.number_of_nodes()
    m = G.number_of_edges()

    if n <= 3:
        print("Skipping graph: The problem states n > 3.")
        print("-" * (len(graph_name) + 24))
        return

    # Decompose into connected components to find c and b
    components = [G.subgraph(c) for c in nx.connected_components(G)]
    c = len(components)
    b = 0
    for comp in components:
        if nx.is_bipartite(comp):
            b += 1
    
    # Calculate k = null(B^T B) = m - n + b
    k = m - n + b
    
    # Calculate lambda_n(G), the largest eigenvalue of the Laplacian
    # Use 'to_numpy_array' to ensure node order is consistent
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]

    print(f"Number of nodes n = {n}")
    print(f"Number of edges m = {m}")
    print(f"Number of connected components c = {c}")
    print(f"Number of bipartite components b = {b}")
    print(f"Calculated k = m - n + b = {m} - {n} + {b} = {k}")
    print(f"Largest Laplacian eigenvalue λ_n(G) = {lambda_n:.4f}")
    
    # Check if the core condition from the problem holds for G
    # We use a small tolerance for floating point comparison
    condition_holds = abs(k - lambda_n / 2) < 1e-9

    if condition_holds:
        print(f"SUCCESS: The condition k = λ_n/2 holds for this graph ({k} = {lambda_n:.4f}/2).")
        
        # Now, check the possible answers for this graph
        print("\nChecking Answer Choices:")
        
        # Choice A: Hard to verify computationally, so we'll do it manually for one case.
        
        # Choice B: The graph has at least k connected components.
        is_B_true = (c >= k)
        print(f"Choice B (c >= k): Is {c} >= {k}? {is_B_true}.")
        if not is_B_true:
            print(">>> This graph falsifies Choice B.")

        # Choice C: The graph is connected.
        is_C_true = (c == 1)
        print(f"Choice C (c == 1): Is {c} == 1? {is_C_true}.")
        if not is_C_true:
            print(">>> This graph falsifies Choice C.")

        # Choice D: The graph has exactly k connected components.
        is_D_true = (c == k)
        print(f"Choice D (c == k): Is {c} == {k}? {is_D_true}.")
        if not is_D_true:
            print(">>> This graph falsifies Choice D.")
    else:
        print(f"INFO: The condition k = λ_n/2 does NOT hold for this graph ({k} != {lambda_n:.4f}/2).")
        
    print("-" * (len(graph_name) + 24))

# Example 1: The complete graph on 4 vertices, K_4
G1 = nx.complete_graph(4)
analyze_graph_property(G1, "Complete graph K_4")

# Example 2: K_4 with a disjoint isolated vertex
G2 = nx.disjoint_union(nx.complete_graph(4), nx.trivial_graph())
analyze_graph_property(G2, "K_4 union an isolated vertex")

# --- Manual check for Choice A using the K_4 graph ---
# For K_4, k=2. The choice is "If you drop 2 edges, there will be at least two nodes with degree <= 1".
print("--- Manual check of Choice A for K_4 ---")
G_test = nx.complete_graph(4)
# Original degrees are [3, 3, 3, 3].
edges_to_remove = [(0, 1), (2, 3)] # Drop two disjoint edges
G_test.remove_edges_from(edges_to_remove)
degrees_after = [d for n, d in G_test.degree()]
nodes_le_1 = sum(1 for d in degrees_after if d <= 1)
print(f"For K_4, k=2. Dropping edges {edges_to_remove}...")
print(f"Degrees become: {degrees_after}")
print(f"Number of nodes with degree <= 1 is {nodes_le_1}.")
if nodes_le_1 < 2:
    print(">>> This case falsifies Choice A.")
print("------------------------------------------")

# Summary of results
print("\n--- Summary ---")
print("Analysis of K_4 showed that the given property can hold for a graph where Choices B and D are false.")
print("Analysis of (K_4 union v) showed that the property can hold for a graph where Choice C is false.")
print("Manual check on K_4 showed a case where Choice A is false.")
print("Since we have found counterexamples for choices A, B, C, and D, none of them can be the correct answer.")
print("Therefore, the statement implies None of the above.")
