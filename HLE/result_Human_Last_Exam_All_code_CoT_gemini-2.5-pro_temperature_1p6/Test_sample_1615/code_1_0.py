def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph
    with 12345 vertices that is not a complete graph.
    """
    # The number of vertices in the graph
    num_vertices = 12345

    # For any graph with n vertices that is not a complete graph,
    # the maximum number of colors required for a proper vertex coloring is n - 1.
    # This is because a graph needs n colors if and only if it is a complete graph (K_n).
    # Since our graph is not complete, its chromatic number is at most n-1.
    # This maximum is achievable for a graph consisting of a K_{n-1} component
    # and a single isolated vertex.
    max_colors = num_vertices - 1

    # Print the equation and the final answer.
    print(f"{num_vertices} - 1 = {max_colors}")

solve_coloring_problem()