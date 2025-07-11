import networkx as nx

def main():
    """
    This script verifies the algebraic properties of graph operations
    to determine the correct classification among the given choices.
    """
    # Define some example simple graphs
    # G1: Path graph on 3 vertices
    G1 = nx.path_graph(3)
    # G2: Path graph on 2 vertices (a single edge)
    G2 = nx.path_graph(2)
    # G3: Complete graph on 3 vertices (a triangle)
    G3 = nx.complete_graph(3)

    print("--- Analyzing the algebraic structure (G, U, x) ---")

    # Property 1: Is multiplication (tensor product) commutative?
    # We check if G1 x G2 is isomorphic to G2 x G1.
    g1_x_g2 = nx.tensor_product(G1, G2)
    g2_x_g1 = nx.tensor_product(G2, G1)
    is_comm = nx.is_isomorphic(g1_x_g2, g2_x_g1)
    print(f"1. Commutativity of Multiplication (x):")
    print(f"   Is G1 x G2 isomorphic to G2 x G1? {is_comm}")
    if is_comm:
        print("   Conclusion: The tensor product is commutative.\n")
    else:
        print("   Conclusion: The tensor product is not commutative.\n")


    # Property 2: Does multiplication (x) distribute over addition (U)?
    # We check if G1 x (G2 U G3) is isomorphic to (G1 x G2) U (G1 x G3).
    # nx.disjoint_union performs the disjoint union operation.
    g2_u_g3 = nx.disjoint_union(G2, G3)
    lhs = nx.tensor_product(G1, g2_u_g3)

    g1_x_g2_dist = nx.tensor_product(G1, G2)
    g1_x_g3_dist = nx.tensor_product(G1, G3)
    rhs = nx.disjoint_union(g1_x_g2_dist, g1_x_g3_dist)
    is_dist = nx.is_isomorphic(lhs, rhs)
    print(f"2. Distributivity of Multiplication (x) over Addition (U):")
    print(f"   LHS graph (G1 x (G2 U G3)) has {lhs.number_of_nodes()} nodes and {lhs.number_of_edges()} edges.")
    print(f"   RHS graph ((G1 x G2) U (G1 x G3)) has {rhs.number_of_nodes()} nodes and {rhs.number_of_edges()} edges.")
    print(f"   Is LHS isomorphic to RHS? {is_dist}")
    if is_dist:
        print("   Conclusion: The tensor product distributes over disjoint union.\n")
    else:
        print("   Conclusion: Distributivity does not hold.\n")

    # Property 3: Does multiplication (x) have an identity element?
    # The only candidate for an identity in simple graphs is K1 (a single vertex).
    K1 = nx.complete_graph(1)
    g1_x_k1 = nx.tensor_product(G1, K1)
    is_id = nx.is_isomorphic(G1, g1_x_k1)
    print(f"3. Multiplicative Identity for Tensor Product (x):")
    print(f"   G1 has {G1.number_of_nodes()} nodes and {G1.number_of_edges()} edges.")
    print(f"   G1 x K1 has {g1_x_k1.number_of_nodes()} nodes and {g1_x_k1.number_of_edges()} edges.")
    print(f"   Is G1 isomorphic to G1 x K1? {is_id}")
    if not is_id:
        print("   Conclusion: There is NO multiplicative identity for the tensor product.\n")
    else:
        print("   Conclusion: A multiplicative identity exists.\n")
        
    # Property 4: Are there additive inverses for addition (U)?
    # A non-empty graph G cannot have an inverse G' such that G U G' is the empty graph.
    print(f"4. Additive Inverses for Disjoint Union (U):")
    print(f"   The number of vertices in G1 U G2 is {G1.number_of_nodes() + G2.number_of_nodes()}.")
    print(f"   This can only be 0 if both graphs are empty.")
    print(f"   Conclusion: Non-empty graphs do NOT have additive inverses. The structure is NOT a ring.\n")

    print("--- Final Conclusion ---")
    print("The structure (G, U, x) forms a commutative semi-ring (assuming the definition does not require a multiplicative identity), but it is not a ring because it lacks additive inverses.")

if __name__ == "__main__":
    main()
<<<D>>>