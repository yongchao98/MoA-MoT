def solve_coloring_problem():
    """
    Calculates the maximum number of colors for a graph with n vertices
    that is not a complete graph.
    """
    # Number of vertices in the graph G
    n = 12345

    # A graph on n vertices has a chromatic number of n if and only if it is the complete graph Kn.
    # Since the graph G is not a complete graph, its chromatic number must be less than n.
    # The maximum possible chromatic number is therefore n - 1.
    # This maximum is achievable, for example, by a graph consisting of a K_{n-1}
    # and one isolated vertex. This graph has n vertices, is not complete, and
    # requires n - 1 colors.
    max_colors = n - 1
    
    # Print the explanation and the final equation.
    print(f"Let n be the number of vertices, n = {n}.")
    print("The maximum number of colors for a graph on n vertices is n, achieved only by the complete graph Kn.")
    print("Since our graph is not a complete graph, its chromatic number must be at most n - 1.")
    print("This maximum is achievable for a graph like K_{n-1} with an additional isolated vertex.")
    print("Thus, the maximum number of colors is calculated as:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()