def solve_graph_coloring_problem():
    """
    This function calculates the maximum number of colors needed to properly
    color a graph with 12345 vertices that is not a complete graph.
    """
    # The number of vertices in the graph G.
    num_vertices = 12345

    # A graph with n vertices requires n colors if and only if it is the complete graph K_n.
    # Since the graph G is not a complete graph, its chromatic number must be less than n.
    # Therefore, the maximum possible chromatic number is n - 1.
    # This maximum is achievable by a graph that consists of a complete subgraph K_{n-1}
    # and one additional isolated vertex. Such a graph has n vertices, is not a complete graph,
    # and its chromatic number is exactly n - 1.

    max_colors = num_vertices - 1

    # Print the explanation and the final equation.
    print(f"The number of vertices is n = {num_vertices}.")
    print("For a graph with n vertices that is not the complete graph, the maximum number of colors needed is n - 1.")
    print("The calculation is:")
    print(f"{num_vertices} - 1 = {max_colors}")

solve_graph_coloring_problem()