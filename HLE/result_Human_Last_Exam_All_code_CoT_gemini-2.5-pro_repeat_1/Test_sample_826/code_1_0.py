import networkx as nx

def tensor_product(g1, g2):
    """A wrapper for networkx.tensor_product that handles empty graphs correctly."""
    if g1.number_of_nodes() == 0 or g2.number_of_nodes() == 0:
        return nx.Graph()
    return nx.tensor_product(g1, g2)

def disjoint_union(g1, g2):
    """A wrapper for networkx.disjoint_union."""
    return nx.disjoint_union(g1, g2)

# --- Define test graphs ---
G1 = nx.path_graph(2) # K2
G2 = nx.Graph()       # K1
G2.add_node(0)
G3 = nx.path_graph(2) # K2


print("Let U be disjoint union and x be tensor product.")
print("We will test the two possible distributive laws.")

# --- Case 1: Test if U distributes over x. Structure (G, x, U) ---
print("\n--- Testing distributivity for (G, x, U) ---")
print("Law to check: G1 U (G2 x G3) == (G1 U G2) x (G1 U G3)")

lhs_1 = disjoint_union(G1, tensor_product(G2, G3))
rhs_1 = tensor_product(disjoint_union(G1, G2), disjoint_union(G1, G3))
dist_holds_1 = nx.is_isomorphic(lhs_1, rhs_1)

print(f"Is the law satisfied? {dist_holds_1}")
if not dist_holds_1:
    print("Counterexample found:")
    print(f"  LHS (G1 U (G2 x G3)) has {lhs_1.number_of_nodes()} nodes and {lhs_1.number_of_edges()} edges.")
    print(f"  RHS ((G1 U G2) x (G1 U G3)) has {rhs_1.number_of_nodes()} nodes and {rhs_1.number_of_edges()} edges.")
    print("Since the graphs are not isomorphic, U does not distribute over x.")
    print("This rules out options C and E.")

# --- Case 2: Test if x distributes over U. Structure (G, U, x) ---
print("\n--- Testing distributivity for (G, U, x) ---")
print("Law to check: G1 x (G2 U G3) == (G1 x G2) U (G1 x G3)")

lhs_2 = tensor_product(G1, disjoint_union(G2, G3))
rhs_2 = disjoint_union(tensor_product(G1, G2), tensor_product(G1, G3))
dist_holds_2 = nx.is_isomorphic(lhs_2, rhs_2)

print(f"Is the law satisfied? {dist_holds_2}")
if dist_holds_2:
    print("This distributive property holds. This is a requirement for options A, B, D.")


print("\n--- Final Analysis ---")
print("1. (G, x, U) is not a semi-ring (rules out C, E).")
print("2. For (G, U, x), the structure is commutative because the tensor product is commutative (rules out B).")
print("3. For (G, U, x), there are no additive inverses for non-empty graphs, so it is not a ring.")
print("4. The choice is between A (not a semi-ring) and D (is a commutative semi-ring, but not a ring).")
print("5. This depends on whether a multiplicative identity is required. Without this requirement, it is a commutative semi-ring. This is a common definition and fits option D perfectly.")
print("\nConclusion: The most appropriate answer is D.")
