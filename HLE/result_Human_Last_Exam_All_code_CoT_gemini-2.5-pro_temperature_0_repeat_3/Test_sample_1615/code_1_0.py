def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with n vertices that is not a complete graph.
    """
    # Number of vertices in the graph G
    num_vertices = 12345

    # For any graph with n vertices, the chromatic number is at most n.
    # The chromatic number is n if and only if the graph is a complete graph (K_n).
    # Since the graph is not a complete graph, its chromatic number is at most n-1.
    # A graph K_n with one edge removed is not complete and has a clique of size n-1,
    # so its chromatic number is exactly n-1.
    # Therefore, the maximum number of colors is n-1.
    max_colors = num_vertices - 1

    # Print the equation and the result
    print(f"The number of vertices is n = {num_vertices}.")
    print("The graph is not a complete graph.")
    print("The maximum number of colors required is n - 1.")
    print(f"Calculation: {num_vertices} - 1 = {max_colors}")
    print(f"So, the maximum number of colours is {max_colors}.")

solve_coloring_problem()