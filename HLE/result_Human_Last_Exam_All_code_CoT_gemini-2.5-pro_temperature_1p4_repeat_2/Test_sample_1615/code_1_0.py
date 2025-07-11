def solve_coloring_problem():
    """
    Calculates the maximum number of colors for a proper vertex coloring of a graph G
    with 12345 vertices, given that G is not a complete graph.
    """
    num_vertices = 12345

    # A fundamental theorem in graph theory states that for a graph with n vertices,
    # its chromatic number (the minimum number of colors for a proper coloring)
    # is equal to n if and only if the graph is a complete graph (K_n),
    # where every vertex is connected to every other vertex.

    # The problem specifies that the graph G is not a complete graph.
    # Therefore, its chromatic number must be less than the total number of vertices.
    max_colors = num_vertices - 1

    print("Let n be the number of vertices in the graph G.")
    print(f"Given n = {num_vertices}")
    print("\nThe problem states that the graph G is not a complete graph.")
    print("\nA key result from graph theory states that the chromatic number of a graph with n vertices, χ(G),")
    print("is equal to n if and only if G is the complete graph K_n.")
    print("\nSince G is not a complete graph, its chromatic number must be strictly less than n.")
    print(f"This means the maximum possible chromatic number is n - 1.")
    print("\nTo confirm this maximum is achievable, consider a graph formed by a complete subgraph K_{n-1} and one isolated vertex.")
    print(f"This graph has {num_vertices} vertices, is not complete, and requires {max_colors} colors to be properly colored.")
    print("\nTherefore, we calculate the maximum number of colors needed as:")
    print(f"{num_vertices} - 1 = {max_colors}")
    print(f"\nThe maximum number of colours that we need to colour the vertices of G properly is {max_colors}.")

solve_coloring_problem()