import networkx as nx

def print_graph_info(G, name):
    """Prints the number of nodes and edges of a graph."""
    nodes = G.number_of_nodes()
    edges = G.number_of_edges()
    print(f"Graph '{name}': Nodes = {nodes}, Edges = {edges}")

# --- Define Graphs for Demonstration ---
# G1 = P3 (Path graph on 3 vertices)
G1 = nx.path_graph(3)
# G2 = K2 (Complete graph on 2 vertices, i.e., a single edge)
G2 = nx.path_graph(2)
# G3 = C3 (Cycle graph on 3 vertices, a triangle)
G3 = nx.cycle_graph(3)
# K0 = Additive identity (empty graph)
K0 = nx.empty_graph(0)
# K1 = Potential multiplicative identity (single vertex graph)
K1 = nx.Graph()
K1.add_node(0)


print("Analysis of the algebraic structure for simple graphs (G)")
print("with disjoint union (U) and tensor product (X).\n")

# --- Analysis of (G, U, X) ---
print("--- Case 1: (G, U, X) where U is addition, X is multiplication ---")

# 1. Check Additive Properties (G, U)
print("\n1. Checking Additive Structure (G, U):")
union_g1_k0 = nx.disjoint_union(G1, K0)
is_identity_add = nx.is_isomorphic(G1, union_g1_k0)
print(f"Does K0 act as additive identity? (P3 U K0 == P3): {is_identity_add}")
print("Disjoint union is also commutative and associative.")
print("Conclusion: (G, U) is a commutative monoid.")

# 2. Check Multiplicative Properties (G, X)
print("\n2. Checking Multiplicative Structure (G, X):")
tensor_g1_g2 = nx.tensor_product(G1, G2)
tensor_g2_g1 = nx.tensor_product(G2, G1)
is_commutative_mul = nx.is_isomorphic(tensor_g1_g2, tensor_g2_g1)
print(f"Is tensor product commutative? (P3 X K2 == K2 X P3): {is_commutative_mul}")

print("Checking for a multiplicative identity...")
tensor_g1_k1 = nx.tensor_product(G1, K1)
is_identity_mul = nx.is_isomorphic(G1, tensor_g1_k1)
print(f"Does K1 act as multiplicative identity? (P3 X K1 == P3): {is_identity_mul}")
print_graph_info(G1, "P3")
print_graph_info(tensor_g1_k1, "P3 X K1")
print("P3 has 2 edges, but P3 X K1 has 0 edges. They are not isomorphic.")
print("Conclusion: (G, X) is a commutative semigroup but NOT a monoid (no identity).")

# 3. Check Distributivity
print("\n3. Checking Distributivity of X over U:")
lhs = nx.tensor_product(G2, nx.disjoint_union(G1, G3))
rhs = nx.disjoint_union(nx.tensor_product(G2, G1), nx.tensor_product(G2, G3))
is_distributive = nx.is_isomorphic(lhs, rhs)
print(f"Does X distribute over U? (K2 X (P3 U C3) == (K2 X P3) U (K2 X C3)): {is_distributive}")

# 4. Check for Ring Property
print("\n4. Checking for Ring Property (Additive Inverses):")
print("For a non-empty graph G, we need an inverse H such that G U H = K0.")
g1_nodes = G1.number_of_nodes()
k0_nodes = K0.number_of_nodes()
print(f"Number of nodes in P3 is {g1_nodes}. Number of nodes in K0 is {k0_nodes}.")
print(f"For any graph H, nodes(P3 U H) = {g1_nodes} + nodes(H), which cannot be {k0_nodes}.")
print("Conclusion: Additive inverses do not exist. The structure is not a ring.")

# 5. Final Conclusion
print("\n--- Final Conclusion ---")
print("The structure (G, U, X) satisfies all properties of a commutative semi-ring")
print("EXCEPT for the existence of a multiplicative identity.")
print("If 'semi-ring' requires a multiplicative identity, the answer is A (not a semi-ring).")
print("If 'semi-ring' does not (a common variant), it is a commutative semi-ring but not a ring.")
print("Given the options, the latter interpretation is the most plausible.")
print("This corresponds to option D.")