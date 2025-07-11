import networkx as nx

def graphs_equal(g1, g2):
    """Checks if two graphs are isomorphic."""
    return nx.is_isomorphic(g1, g2)

# Define three simple graphs for the demonstration.
# We use path_graph(n) which creates a simple path of n vertices.
# G1 = P_2 (isomorphic to K_2)
# G2 = P_2 (isomorphic to K_2)
# G3 = P_1 (isomorphic to K_1)
G1 = nx.path_graph(2)
G2 = nx.path_graph(2)
G3 = nx.path_graph(1)

# It's good practice to use disjoint node sets for verification, though networkx handles it.
G1 = nx.relabel_nodes(G1, {0: 'g1_0', 1: 'g1_1'})
G2 = nx.relabel_nodes(G2, {0: 'g2_0', 1: 'g2_1'})
G3 = nx.relabel_nodes(G3, {0: 'g3_0'})

print("We will verify the distributivity law for graph operations:")
print("G1 tensor (G2 union G3) is isomorphic to (G1 tensor G2) union (G1 tensor G3)")
print(f"Let G1 be P_2 (2 vertices, 1 edge)")
print(f"Let G2 be P_2 (2 vertices, 1 edge)")
print(f"Let G3 be P_1 (1 vertex, 0 edges)\n")

# --- Left Hand Side (LHS) ---
# First, compute the disjoint union of G2 and G3
G2_union_G3 = nx.disjoint_union(G2, G3)
# Then, compute the tensor product with G1
LHS = nx.tensor_product(G1, G2_union_G3)

# --- Right Hand Side (RHS) ---
# First, compute the tensor products
G1_tensor_G2 = nx.tensor_product(G1, G2)
G1_tensor_G3 = nx.tensor_product(G1, G3)
# Then, compute the disjoint union of the results
RHS = nx.disjoint_union(G1_tensor_G2, G1_tensor_G3)

print("--- Calculating LHS: G1 tensor (G2 union G3) ---")
print(f"Number of vertices in G2 union G3: {G2_union_G3.number_of_vertices()}")
print(f"Number of edges in G2 union G3: {G2_union_G3.number_of_edges()}")
print(f"Resulting LHS graph has {LHS.number_of_vertices()} vertices and {LHS.number_of_edges()} edges.\n")

print("--- Calculating RHS: (G1 tensor G2) union (G1 tensor G3) ---")
print(f"Number of vertices in G1 tensor G2: {G1_tensor_G2.number_of_vertices()}")
print(f"Number of edges in G1 tensor G2: {G1_tensor_G2.number_of_edges()}")
print(f"Number of vertices in G1 tensor G3: {G1_tensor_G3.number_of_vertices()}")
print(f"Number of edges in G1 tensor G3: {G1_tensor_G3.number_of_edges()}")
print(f"Resulting RHS graph has {RHS.number_of_vertices()} vertices and {RHS.number_of_edges()} edges.\n")

# --- Final Verification ---
are_isomorphic = graphs_equal(LHS, RHS)
print(f"Verification: Are the LHS and RHS graphs isomorphic? {are_isomorphic}")
if are_isomorphic:
    print("The distributivity property holds true for this example.")
