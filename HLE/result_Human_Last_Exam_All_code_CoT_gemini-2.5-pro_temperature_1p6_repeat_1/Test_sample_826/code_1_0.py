import networkx as nx

def demonstrate_distributivity():
    """
    This function demonstrates the distributive property of the tensor product
    over the disjoint union of graphs, i.e., G1 @ (G2 U G3) ~= (G1 @ G2) U (G1 @ G3),
    where @ is tensor product and U is disjoint union.
    We check this property using sample graphs from the networkx library.
    This property is a key requirement for (G, U, @) to be a semi-ring.
    """
    # 1. Define three simple graphs
    # G1: Path graph on 3 vertices (P3)
    # G2: Complete graph on 2 vertices (K2, a single edge)
    # G3: A single isolated vertex (K1)
    G1 = nx.path_graph(3)
    G2 = nx.complete_graph(2)
    G3 = nx.complete_graph(1)

    print("Demonstrating the distributive property for simple graphs.")
    print("Let '+' be disjoint union (U) and '*' be tensor product (@).")
    print("We will verify if G1 * (G2 + G3) is isomorphic to (G1 * G2) + (G1 * G3).")
    print(f"G1: {G1.nodes}, Edges: {G1.edges}")
    print(f"G2: {G2.nodes}, Edges: {G2.edges}")
    print(f"G3: {G3.nodes}, Edges: {G3.edges}")
    print("-" * 20)

    # 2. Compute the left-hand side: G1 * (G2 + G3)
    G2_union_G3 = nx.disjoint_union(G2, G3)
    LHS = nx.tensor_product(G1, G2_union_G3)
    print("Left Hand Side Graph (LHS): G1 @ (G2 U G3)")
    print(f"Number of nodes: {LHS.number_of_nodes()}")
    print(f"Number of edges: {LHS.number_of_edges()}")

    # 3. Compute the right-hand side: (G1 * G2) + (G1 * G3)
    G1_tensor_G2 = nx.tensor_product(G1, G2)
    G1_tensor_G3 = nx.tensor_product(G1, G3)
    RHS = nx.disjoint_union(G1_tensor_G2, G1_tensor_G3)
    print("\nRight Hand Side Graph (RHS): (G1 @ G2) U (G1 @ G3)")
    print(f"Number of nodes: {RHS.number_of_nodes()}")
    print(f"Number of edges: {RHS.number_of_edges()}")
    print("-" * 20)

    # 4. Check for isomorphism
    are_isomorphic = nx.is_isomorphic(LHS, RHS)

    print(f"Are LHS and RHS graphs isomorphic? {are_isomorphic}")

    if are_isomorphic:
        print("\nThe distributivity property holds for our example graphs.")
        print("This supports the conclusion that (G, U, @) satisfies the semi-ring axioms.")
    else:
        print("\nThe distributivity property does not hold for our example.")

if __name__ == '__main__':
    # Running the demonstration function
    # Note: You may need to install networkx (`pip install networkx`)
    try:
        demonstrate_distributivity()
    except ImportError:
        print("Please install the networkx library to run this demonstration (`pip install networkx`).")
