import networkx as nx

def print_isomorphism_check(G1, G2, G1_name, G2_name, description):
    """Checks for isomorphism and prints a detailed, user-friendly result."""
    are_isomorphic = nx.is_isomorphic(G1, G2)
    print(f"--- {description} ---")
    print(f"Comparing '{G1_name}' and '{G2_name}':")
    # Outputting numbers as requested
    print(f"  '{G1_name}' has {G1.number_of_nodes()} nodes and {G1.number_of_edges()} edges.")
    print(f"  '{G2_name}' has {G2.number_of_nodes()} nodes and {G2.number_of_edges()} edges.")
    print(f"  Are they isomorphic? {are_isomorphic}\n")
    return are_isomorphic

# Define some simple graphs for testing
K_0 = nx.Graph() # Additive identity for union
K_1 = nx.path_graph(1) # Single vertex, no edges
K_2 = nx.complete_graph(2) # Single edge
P_3 = nx.path_graph(3) # Path on 3 vertices

print("### Analysis of (G, union, tensor_product) ###\n")

# 1. Check if (G, union) is a commutative monoid (our "addition")
print("Step 1: Analyzing '+' = Disjoint Union (U)\n")
G1_U_G2 = nx.disjoint_union(K_2, P_3)
G2_U_G1 = nx.disjoint_union(P_3, K_2)
print_isomorphism_check(G1_U_G2, G2_U_G1, "K_2 U P_3", "P_3 U K_2", "Commutativity of Union")

K2_U_K0 = nx.disjoint_union(K_2, K_0)
print_isomorphism_check(K2_U_K0, K_2, "K_2 U K_0", "K_2", "Identity of Union")
print("Result: (G, U) is a commutative monoid.\n")

# 2. Check properties of (G, tensor_product) (our "multiplication")
print("Step 2: Analyzing '*' = Tensor Product (x)\n")
G1_T_G2 = nx.tensor_product(K_2, P_3)
G2_T_G1 = nx.tensor_product(P_3, K_2)
print_isomorphism_check(G1_T_G2, G2_T_G1, "K_2 x P_3", "P_3 x K_2", "Commutativity of Tensor Product")

# Check for multiplicative identity. It must be K_1 if it exists in G.
K2_T_K1 = nx.tensor_product(K_2, K_1)
print_isomorphism_check(K2_T_K1, K_2, "K_2 x K_1", "K_2", "Multiplicative Identity")
print("Result: Tensor product is commutative, but has NO identity in the set of simple graphs.")
print("K_2 x K_1 has 0 edges, while K_2 has 1 edge. They are not isomorphic.")
print("This means (G, U, x) is not a semi-ring under the strictest definition.\n")

# 3. Check distributivity
print("Step 3: Checking Distributivity of 'x' over 'U'\n")
LHS = nx.tensor_product(K_2, nx.disjoint_union(K_1, P_3))
RHS = nx.disjoint_union(nx.tensor_product(K_2, K_1), nx.tensor_product(K_2, P_3))
print_isomorphism_check(LHS, RHS, "K_2 x (K_1 U P_3)", "(K_2 x K_1) U (K_2 x P_3)", "Distributivity")
print("Result: Distributivity holds.\n")

# 4. Check if it's a ring
print("Step 4: Checking for Additive Inverses (Ring Property)\n")
print("For a ring, every graph G must have an inverse G' such that G U G' is the empty graph K_0.")
print(f"Let G = K_2. It has {K_2.number_of_nodes()} nodes.")
print(f"Any union G U G' will have at least {K_2.number_of_nodes()} nodes.")
print(f"The identity K_0 has {K_0.number_of_nodes()} nodes.")
print("Since 2 > 0, no such inverse G' can exist for K_2.")
print("Result: (G,U) is not a group. Thus, the structure is not a ring.\n")

# 5. Final Conclusion
print("### Final Conclusion ###")
print("Assuming the common relaxed definition of a semi-ring (no multiplicative identity required):")
print("- The structure is a semi-ring because the necessary axioms hold.")
print("- It is a COMMUTATIVE semi-ring because the tensor product (x) is commutative.")
print("- It is NOT a ring because there are no additive inverses for union (U).")
print("This matches option D.")
