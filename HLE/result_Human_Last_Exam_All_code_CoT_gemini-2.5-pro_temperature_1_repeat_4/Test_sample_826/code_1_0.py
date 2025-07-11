import networkx as nx

def main():
    """
    This script checks the distributivity law for the two potential semi-ring structures.
    It demonstrates that Union does NOT distribute over Tensor Product, eliminating
    (G, tensor, union) as a candidate for a semi-ring.
    """
    # Create the graphs K_1 and K_2
    G1 = nx.Graph()
    G1.add_node(0)  # K_1

    G2 = nx.path_graph(2)  # K_2
    G3 = nx.path_graph(2)  # K_2, labelled differently for clarity
    nx.relabel_nodes(G3, {0: 'a', 1: 'b'}, copy=False)

    print("We test the distributivity law G1 U (G2 x G3) == (G1 U G2) x (G1 U G3)")
    print("where U is disjoint union and x is tensor product.")
    print("Let G1 = K_1 (1 vertex, 0 edges) and G2 = G3 = K_2 (2 vertices, 1 edge).\n")

    # --- LHS: G1 U (G2 x G3) ---
    print("--- Left Hand Side: G1 U (G2 x G3) ---")
    # Calculate G2 x G3
    G2_x_G3 = nx.tensor_product(G2, G3)
    num_nodes_G2_x_G3 = G2_x_G3.number_of_nodes()
    num_edges_G2_x_G3 = G2_x_G3.number_of_edges()
    print(f"The graph G2 x G3 has {num_nodes_G2_x_G3} vertices and {num_edges_G2_x_G3} edges.")

    # Calculate G1 U (G2 x G3)
    LHS = nx.disjoint_union(G1, G2_x_G3)
    num_nodes_LHS = LHS.number_of_nodes()
    num_edges_LHS = LHS.number_of_edges()

    print(f"The graph G1 U (G2 x G3) has {G1.number_of_nodes()} + {num_nodes_G2_x_G3} = {num_nodes_LHS} vertices.")
    print(f"The graph G1 U (G2 x G3) has {G1.number_of_edges()} + {num_edges_G2_x_G3} = {num_edges_LHS} edges.\n")


    # --- RHS: (G1 U G2) x (G1 U G3) ---
    print("--- Right Hand Side: (G1 U G2) x (G1 U G3) ---")
    # Calculate G1 U G2
    G1_U_G2 = nx.disjoint_union(G1, G2)
    num_nodes_G1_U_G2 = G1_U_G2.number_of_nodes()
    num_edges_G1_U_G2 = G1_U_G2.number_of_edges()
    print(f"The graph G1 U G2 has {num_nodes_G1_U_G2} vertices and {num_edges_G1_U_G2} edges.")

    # In our case, G1 U G3 is isomorphic to G1 U G2
    G1_U_G3 = nx.disjoint_union(G1, G3)
    
    # Calculate (G1 U G2) x (G1 U G3)
    RHS = nx.tensor_product(G1_U_G2, G1_U_G3)
    num_nodes_RHS = RHS.number_of_nodes()
    num_edges_RHS = RHS.number_of_edges()

    print(f"The graph (G1 U G2) x (G1 U G3) has {num_nodes_G1_U_G2} * {num_nodes_G1_U_G2} = {num_nodes_RHS} vertices.")
    # The number of edges of a tensor product GxH is not a simple product of edge counts.
    # We let networkx compute it for us.
    print(f"The graph (G1 U G2) x (G1 U G3) has {num_edges_RHS} edges.\n")


    # --- Comparison ---
    print("--- Comparison ---")
    print(f"LHS has {num_nodes_LHS} vertices and {num_edges_LHS} edges.")
    print(f"RHS has {num_nodes_RHS} vertices and {num_edges_RHS} edges.")
    if num_nodes_LHS != num_nodes_RHS or num_edges_LHS != num_edges_RHS:
        print("Since the number of vertices is different, the graphs are not isomorphic.")
        print("Therefore, U does not distribute over x.")
        print("This rules out (G, x, U) as a semi-ring.")
    else:
        print("The number of vertices and edges are the same. This test is inconclusive.")

if __name__ == "__main__":
    main()
