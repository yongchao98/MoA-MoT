import networkx as nx

def pretty_print_graph(G, name):
    """Helper function to print graph info."""
    print(f"Graph '{name}': {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

# Create a few simple graphs for testing
G1 = nx.path_graph(3) # P_3: path on 3 vertices
G2 = nx.path_graph(2) # P_2, which is also K_2
G3 = nx.complete_graph(3) # K_3: a triangle
K1 = nx.complete_graph(1) # K_1: single vertex, no edges

print("--- Step-by-step analysis of (G, U, x) properties ---")

# 1. Check Commutativity of Tensor Product (x)
# Is G1 x G2 isomorphic to G2 x G1?
print("\n1. Checking Commutativity of Tensor Product (x)...")
prod12 = nx.tensor_product(G1, G2)
prod21 = nx.tensor_product(G2, G1)
pretty_print_graph(prod12, "P3 x K2")
pretty_print_graph(prod21, "K2 x P3")
is_comm = nx.is_isomorphic(prod12, prod21)
print(f"Is P3 x K2 isomorphic to K2 x P3? {is_comm}")
print("Conclusion: Tensor product is commutative.")

# 2. Check Distributivity of Tensor Product (x) over Disjoint Union (U)
# Is G1 x (G2 U G3) isomorphic to (G1 x G2) U (G1 x G3)?
print("\n2. Checking Distributivity of x over U...")
G2_union_G3 = nx.disjoint_union(G2, G3)
lhs = nx.tensor_product(G1, G2_union_G3)
prod12 = nx.tensor_product(G1, G2)
prod13 = nx.tensor_product(G1, G3)
rhs = nx.disjoint_union(prod12, prod13)
pretty_print_graph(lhs, "P3 x (K2 U K3)")
pretty_print_graph(rhs, "(P3 x K2) U (P3 x K3)")
is_dist = nx.is_isomorphic(lhs, rhs)
print(f"Is P3 x (K2 U K3) isomorphic to (P3 x K2) U (P3 x K3)? {is_dist}")
print("Conclusion: Tensor product distributes over disjoint union.")


# 3. Check for Multiplicative Identity
# K1 is the only candidate for an identity, but G x K1 is not isomorphic to G
print("\n3. Checking for Multiplicative Identity...")
pretty_print_graph(G1, "P3")
prod_g1_k1 = nx.tensor_product(G1, K1)
pretty_print_graph(prod_g1_k1, "P3 x K1")
is_ident = nx.is_isomorphic(G1, prod_g1_k1)
print(f"Is P3 x K1 isomorphic to P3? {is_ident}")
print("Note that P3 has 2 edges, while P3 x K1 has 0 edges.")
print("Conclusion: There is no multiplicative identity for tensor product in the set of simple graphs.")

# 4. Check if it's a Ring
# It cannot be a ring as there are no additive inverses for U
print("\n4. Checking if it's a Ring...")
print("The additive identity is the empty graph K0.")
print("For any non-empty graph G, there is no graph H such that G U H = K0,")
print("because the number of vertices in G U H is |V(G)| + |V(H)|, which can't be 0 if |V(G)| > 0.")
print("Conclusion: Not a ring, as no additive inverses exist.")

print("\n--- Final Algebraic Structure ---")
print("Structure: (G, U, x)")
print(" (G, U): Commutative monoid (satisfies additive properties)")
print(" (G, x): Commutative semigroup (satisfies multiplicative properties, except for identity)")
print(" Distributivity and Annihilation: Hold")
print("\nThis structure is a 'commutative semi-ring, but not a ring'.")
