import numpy as np
import networkx as nx

def check_graph_property(G):
    """
    Checks if a graph G satisfies the property:
    null(B^T B) = lambda_n(G) / 2
    which is equivalent to:
    m - n + k_b = lambda_n(L) / 2
    """
    if G.number_of_nodes() == 0:
        print("Graph is empty.")
        return

    # 1. Get n, m, k
    n = G.number_of_nodes()
    m = G.number_of_edges()
    k = nx.number_connected_components(G)
    print(f"Graph properties: n={n}, m={m}, k={k}")

    # 2. Get k_b (number of bipartite components)
    k_b = 0
    components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    for i, component in enumerate(components):
        if nx.is_bipartite(component):
            k_b += 1
            # print(f"  - Component {i+1} is bipartite.")
        # else:
            # print(f"  - Component {i+1} is not bipartite.")
    print(f"Number of bipartite components (k_b): {k_b}")

    # 3. Calculate LHS: m - n + k_b
    lhs = m - n + k_b
    print(f"LHS of the equation (m - n + k_b) = {m} - {n} + {k_b} = {lhs}")

    # 4. Get lambda_n(G)
    # The graph can be disconnected, networkx.laplacian_spectrum handles this correctly.
    try:
        eigenvalues = sorted(nx.laplacian_spectrum(G))
        lambda_n = eigenvalues[-1]
    except nx.NetworkXError:
        # Handle case with no edges
        if m == 0:
            lambda_n = 0
        else:
            raise
            
    print(f"Largest eigenvalue of Laplacian (lambda_n): {lambda_n:.4f}")

    # 5. Calculate RHS: lambda_n / 2
    rhs = lambda_n / 2
    print(f"RHS of the equation (lambda_n / 2) = {lambda_n:.4f} / 2 = {rhs:.4f}")
    
    # 6. Check if the property holds
    if np.isclose(lhs, rhs):
        print("\nConclusion: The property HOLDS for this graph.")
    else:
        print("\nConclusion: The property DOES NOT HOLD for this graph.")
    
    # 7. Check options B, C, D
    print("\n--- Checking Answer Choices ---")
    # Option C
    is_connected = (k == 1)
    print(f"Option C (Graph is connected): {is_connected}")
    # Option B
    holds_b = (k >= rhs)
    print(f"Option B (k >= lambda_n/2 -> {k} >= {rhs:.4f}): {holds_b}")
    # Option D
    holds_d = np.isclose(k, rhs)
    print(f"Option D (k == lambda_n/2 -> {k} == {rhs:.4f}): {holds_d}")


# --- Test Case 1: K4 ---
print("--- Analyzing graph K4 ---")
G1 = nx.complete_graph(4)
check_graph_property(G1)
print("="*30)

# --- Test Case 2: K4 U K1 ---
print("\n--- Analyzing graph K4 U K1 ---")
G2 = nx.complete_graph(4)
G2.add_node(4) # networkx nodes start at 0, so this is the 5th node
check_graph_property(G2)
print("="*30)
