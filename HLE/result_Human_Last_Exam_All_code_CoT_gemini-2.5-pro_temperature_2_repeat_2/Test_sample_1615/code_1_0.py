def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with n vertices, given that it is not a complete graph.
    """
    # Number of vertices in the graph G
    num_vertices = 12345

    # For any graph with n vertices, the maximum number of colors needed is n.
    # This maximum is only achieved by the complete graph K_n.
    # The problem states that our graph G is not a complete graph.
    # Therefore, the number of colors must be strictly less than n.
    # The maximum possible number of colors is n - 1.
    max_colors = num_vertices - 1

    # We can prove this is achievable by constructing a graph.
    # Consider a graph formed by a complete subgraph on (n-1) vertices (a K_{n-1})
    # and one additional isolated vertex. This graph has n vertices, is not complete,
    # and requires (n-1) colors.

    print(f"The number of vertices in the graph is n = {num_vertices}.")
    print("A graph on n vertices requires n colors if and only if it is the complete graph K_n.")
    print("Since the graph G is not the complete graph, it must require fewer than n colors.")
    print("The maximum possible number of colors for a non-complete graph on n vertices is n - 1.")
    print(f"So, the maximum number of colors needed is {num_vertices} - 1 = {max_colors}.")

solve_coloring_problem()