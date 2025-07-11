def solve():
    """
    Computes the first l2-Betti number of the fundamental group G.

    The problem simplifies to computing the first Betti number of the underlying graph,
    which is the line graph of the Petersen graph.
    """

    # Properties of the Petersen graph (P)
    petersen_vertices = 10
    petersen_edges = 15
    petersen_degree = 3  # The Petersen graph is 3-regular

    # Properties of the Line Graph of the Petersen graph (L(P))
    # The number of vertices in L(P) is the number of edges in P.
    line_graph_vertices = petersen_edges

    # The number of edges in L(P) is the sum of combinations of edges meeting at each vertex of P.
    # For a k-regular graph with n vertices, this is n * (k choose 2).
    line_graph_edges = petersen_vertices * (petersen_degree * (petersen_degree - 1) // 2)

    # The first Betti number of a connected graph is given by |E| - |V| + 1.
    # The line graph of a connected graph (like the Petersen graph) is connected.
    betti_number_1 = line_graph_edges - line_graph_vertices + 1

    print("Step 1: The problem reduces to computing the first Betti number of the line graph of the Petersen graph.")
    print(f"Let L(P) be the line graph of the Petersen graph.")
    print(f"The first l2-Betti number is equal to b1(L(P)).\n")

    print("Step 2: Calculate the number of vertices and edges of L(P).")
    print(f"Number of vertices in L(P) = Number of edges in Petersen graph = {line_graph_vertices}")
    print(f"Number of edges in L(P) = {line_graph_edges}\n")

    print("Step 3: Calculate the first Betti number of L(P).")
    print(f"b1(L(P)) = |E(L(P))| - |V(L(P))| + 1")
    print(f"b1(L(P)) = {line_graph_edges} - {line_graph_vertices} + 1 = {betti_number_1}\n")

    print(f"The first l2-Betti number of the fundamental group of G is {betti_number_1}.")


solve()