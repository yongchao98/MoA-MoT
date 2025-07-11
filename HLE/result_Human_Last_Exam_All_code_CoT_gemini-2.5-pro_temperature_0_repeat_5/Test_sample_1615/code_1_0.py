def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with n vertices that is not a complete graph.
    """
    # Number of vertices in the graph G
    n = 12345

    # If a graph on n vertices is not complete, it means there is at least one pair of
    # non-adjacent vertices. These two vertices can be given the same color.
    # The remaining n-2 vertices can be colored with n-2 other distinct colors in the worst case.
    # This gives an upper bound of n-1 colors.
    # To show this is the maximum, consider a complete graph Kn with one edge removed.
    # This graph is not complete, but it still contains a clique of size n-1,
    # which requires n-1 colors.
    # Therefore, the maximum number of colors is n - 1.
    max_colors = n - 1

    # Print the equation as requested
    print(f"The maximum number of colors is given by the equation:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()