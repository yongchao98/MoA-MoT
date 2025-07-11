import networkx as nx

def check_isomorphism(G1, G2):
    """Checks for isomorphism between two graphs."""
    return nx.is_isomorphic(G1, G2)

def main():
    """
    Analyzes the algebraic properties of graph operations to determine the correct structure.
    """
    # Define some simple test graphs
    G0 = nx.Graph() # Empty graph, potential additive identity for Union
    G1 = nx.path_graph(1)  # K1, single vertex
    K2 = nx.path_graph(2)  # K2, single edge
    P3 = nx.path_graph(3)  # Path on 3 vertices

    print("--- Testing Structure 1: (G, Union, Tensor Product) ---")
    print("Here, Addition (+) is Disjoint Union, Multiplication (*) is Tensor Product.")
    
    # Test (G, Union) as a commutative monoid
    print("\n1. Properties of (G, +) = (G, Union):")
    # Commutativity: G1 + G2 == G2 + G1
    add_comm = check_isomorphism(nx.disjoint_union(K2, P3), nx.disjoint_union(P3, K2))
    print(f"  - Commutativity: {add_comm}")
    # Associativity: (G1 + G2) + G3 == G1 + (G2 + G3)
    add_assoc_lhs = nx.disjoint_union(nx.disjoint_union(G1, K2), P3)
    add_assoc_rhs = nx.disjoint_union(G1, nx.disjoint_union(K2, P3))
    add_assoc = check_isomorphism(add_assoc_lhs, add_assoc_rhs)
    print(f"  - Associativity: {add_assoc}")
    # Identity: G + 0 == G
    add_id = check_isomorphism(nx.disjoint_union(P3, G0), P3)
    print(f"  - Has Additive Identity (Empty Graph): {add_id}")

    # Test (G, Tensor) properties
    print("\n2. Properties of (G, *) = (G, Tensor Product):")
    # Commutativity: G1 * G2 == G2 * G1
    mul_comm = check_isomorphism(nx.tensor_product(K2, P3), nx.tensor_product(P3, K2))
    print(f"  - Commutativity: {mul_comm}")
    # Associativity: (G1 * G2) * G3 == G1 * (G2 * G3)
    mul_assoc_lhs = nx.tensor_product(nx.tensor_product(G1, K2), P3)
    mul_assoc_rhs = nx.tensor_product(G1, nx.tensor_product(K2, P3))
    mul_assoc = check_isomorphism(mul_assoc_lhs, mul_assoc_rhs)
    print(f"  - Associativity: {mul_assoc}")
    # Identity: Is there a '1' s.t. G * 1 == G? Try K1.
    mul_id_test = nx.tensor_product(K2, G1)
    mul_id = check_isomorphism(mul_id_test, K2)
    print(f"  - Has Multiplicative Identity: {mul_id} (Test with K2*K1 shows it fails)")
    
    # Test Distributivity: G1 * (G2 + G3) == (G1 * G2) + (G1 * G3)
    print("\n3. Testing Distributivity of Tensor over Union:")
    dist_lhs = nx.tensor_product(K2, nx.disjoint_union(P3, G1))
    dist_rhs = nx.disjoint_union(nx.tensor_product(K2, P3), nx.tensor_product(K2, G1))
    dist = check_isomorphism(dist_lhs, dist_rhs)
    print(f"  - Distributivity: {dist}")
    
    # Test Ring property: existence of additive inverse
    print("\n4. Testing for Additive Inverses (Ring property):")
    # Can we find InvG such that P3 + InvG == G0?
    # P3 has 3 nodes. Union adds nodes. Cannot result in 0 nodes.
    print("  - Has Additive Inverses: False (Cannot return to the empty graph via union)")

    print("\nSummary for (G, Union, Tensor):")
    print("  - It is a commutative structure where multiplication distributes over addition.")
    print("  - It is not a ring.")
    print("  - It lacks a multiplicative identity, so it is not a semi-ring under the strictest definition.")
    print("  - It fits the description of a commutative semi-ring (if mult. identity is not required) but not a ring.")
    print("  - This matches Option D.")

    print("\n\n--- Testing Structure 2: (G, Tensor Product, Union) ---")
    print("Here, Addition (+) is Tensor Product, Multiplication (*) is Disjoint Union.")
    # Test (G, Tensor) as a commutative monoid
    print("\n1. Properties of (G, +) = (G, Tensor Product):")
    # Identity: Is there a '0' s.t. G + 0 == G? (Same as finding mult. identity above)
    add_id_2 = check_isomorphism(nx.tensor_product(K2, G1), K2)
    print(f"  - Has Additive Identity: {add_id_2} (Fails, so not a monoid)")
    
    # Test Distributivity: G1 * (G2 + G3) == (G1 * G2) + (G1 * G3)
    print("\n2. Testing Distributivity of Union over Tensor:")
    dist2_lhs = nx.disjoint_union(G1, nx.tensor_product(K2, K2))
    dist2_rhs = nx.tensor_product(nx.disjoint_union(G1, K2), nx.disjoint_union(G1, K2))
    dist2 = check_isomorphism(dist2_lhs, dist2_rhs)
    print(f"  - Distributivity: {dist2} (Fails based on vertex count: {dist2_lhs.number_of_nodes()} vs {dist2_rhs.number_of_nodes()})")

    print("\nSummary for (G, Tensor, Union):")
    print("  - This structure is not a semi-ring because (G, Tensor Product) is not a monoid, and distributivity fails.")
    print("  - This eliminates options C and E.")
    
if __name__ == "__main__":
    main()