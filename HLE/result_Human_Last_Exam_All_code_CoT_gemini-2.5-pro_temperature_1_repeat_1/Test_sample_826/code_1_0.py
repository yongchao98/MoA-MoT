import networkx as nx

# Create some simple graphs for testing
# P2 is the path graph with 2 vertices (equivalent to K2)
P2 = nx.path_graph(2) 
# P3 is the path graph with 3 vertices
P3 = nx.path_graph(3)
# K1 is the graph with a single vertex
K1 = nx.complete_graph(1)

# --- Analysis of (G, U, O) ---
# O = Tensor Product (@), U = Disjoint Union (+)

print("Analyzing the structure (G, disjoint_union, tensor_product):")

# Check 1: Commutativity of tensor product (multiplication)
# G1 @ G2 is isomorphic to G2 @ G1?
G1 = P2
G2 = P3
is_commutative = nx.is_isomorphic(nx.tensor_product(G1, G2), nx.tensor_product(G2, G1))
print(f"1. Is tensor product commutative? (e.g., P2 @ P3 vs P3 @ P2): {is_commutative}")

# Check 2: Distributivity of tensor product over disjoint union
# G1 @ (G2 + G3) is isomorphic to (G1 @ G2) + (G1 @ G3)?
G3 = K1
# The '+' operator in networkx is not disjoint union, nx.disjoint_union must be used.
LHS_dist = nx.tensor_product(G1, nx.disjoint_union(G2, G3))
RHS_dist = nx.disjoint_union(nx.tensor_product(G1, G2), nx.tensor_product(G1, G3))
is_distributive = nx.is_isomorphic(LHS_dist, RHS_dist)
print(f"2. Does tensor product distribute over disjoint union?: {is_distributive}")

# Check 3: Existence of multiplicative identity
# Is K1 the identity for tensor product? Check if P2 @ K1 is isomorphic to P2.
has_identity = nx.is_isomorphic(nx.tensor_product(P2, K1), P2)
print(f"3. Is K1 the multiplicative identity? (P2 @ K1 == P2): {has_identity}")
p2_edges = P2.number_of_edges()
p2_k1_prod_edges = nx.tensor_product(P2, K1).number_of_edges()
print(f"   - Edges in P2: {p2_edges}, Edges in P2 @ K1: {p2_k1_prod_edges}. They are not isomorphic.")

# Check 4: Existence of additive inverses
# Is there a G' such that G + G' = K0 (empty graph)?
# This is impossible by counting vertices for any non-empty G.
is_a_ring = False
print(f"4. Does it form a ring (are there additive inverses)?: {is_a_ring}")

# --- Analysis of (G, O, U) ---
# O = Tensor Product, U = Disjoint Union
print("\nAnalyzing the structure (G, tensor_product, disjoint_union):")

# Check 5: Distributivity of disjoint union over tensor product
# G1 + (G2 @ G3) is isomorphic to (G1 + G2) @ (G1 + G3)?
LHS_dist2 = nx.disjoint_union(G1, nx.tensor_product(G2, G3))
RHS_dist2 = nx.tensor_product(nx.disjoint_union(G1, G2), nx.disjoint_union(G1, G3))
is_distributive2 = nx.is_isomorphic(LHS_dist2, RHS_dist2)
print(f"5. Does disjoint union distribute over tensor product?: {is_distributive2}")
lhs_nodes = LHS_dist2.number_of_nodes()
rhs_nodes = RHS_dist2.number_of_nodes()
print(f"   - Vertices in LHS graph: {lhs_nodes}")
print(f"   - Vertices in RHS graph: {rhs_nodes}. They are not isomorphic.")
