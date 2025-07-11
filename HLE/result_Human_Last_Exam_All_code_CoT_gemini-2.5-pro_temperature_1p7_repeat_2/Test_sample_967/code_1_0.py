def solve():
    """
    Computes the first l2-Betti number of the fundamental group of G.

    The problem simplifies to computing the first Betti number of the
    line graph of the Petersen graph.
    """

    # Properties of the Petersen graph (P)
    petersen_vertices = 10
    petersen_edges = 15
    # The Petersen graph is a 3-regular graph
    petersen_degree = 3

    # Properties of the line graph of the Petersen graph (L(P))
    # Number of vertices in L(P) is the number of edges in P
    line_graph_vertices = petersen_edges

    # Number of edges in L(P) is sum of (d choose 2) for each vertex in P,
    # where d is the degree of the vertex. Since P is 3-regular:
    line_graph_edges = petersen_vertices * (petersen_degree * (petersen_degree - 1) // 2)

    # The line graph of a connected graph is connected.
    # The first Betti number of a connected graph Gamma is b1 = |E| - |V| + 1
    betti_number = line_graph_edges - line_graph_vertices + 1

    print("Step 1: The problem reduces to computing the first Betti number of the underlying graph,")
    print("which is the line graph of the Petersen graph, L(P).")
    print("b_1(L(P)) = |E(L(P))| - |V(L(P))| + 1\n")

    print(f"Step 2: The number of vertices in L(P) is the number of edges in the Petersen graph.")
    print(f"|V(L(P))| = {line_graph_vertices}\n")

    print(f"Step 3: The number of edges in L(P) is calculated based on the degrees of vertices in the Petersen graph.")
    print(f"The Petersen graph has {petersen_vertices} vertices, each of degree {petersen_degree}.")
    print(f"|E(L(P))| = {petersen_vertices} * C({petersen_degree}, 2) = {petersen_vertices} * {petersen_degree * (petersen_degree - 1) // 2} = {line_graph_edges}\n")

    print("Step 4: Compute the first Betti number.")
    print(f"b_1(L(P)) = |E(L(P))| - |V(L(P))| + 1 = {line_graph_edges} - {line_graph_vertices} + 1 = {betti_number}")
    print("\nThe first l2-betti number is therefore equal to this value.")
    print(f"Final answer: {betti_number}")

solve()