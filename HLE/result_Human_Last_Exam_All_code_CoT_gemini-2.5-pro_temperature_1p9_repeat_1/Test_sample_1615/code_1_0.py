def solve_coloring_problem():
    """
    Calculates the maximum number of colors for a non-complete graph with n vertices.
    """
    # The number of vertices in the graph G
    n = 12345

    # According to graph theory, for any graph G with n vertices that is not a
    # complete graph, the chromatic number chi(G) is at most n-1.
    # This bound is tight, as a graph consisting of a complete graph on n-1
    # vertices and one isolated vertex requires n-1 colors.
    #
    # So, the maximum number of colors is n - 1.

    max_colors = n - 1

    print(f"Let n be the number of vertices, n = {n}.")
    print("The graph is not a complete graph.")
    print("The maximum number of colors required for a proper vertex coloring of such a graph is n - 1.")
    print(f"Calculation: {n} - 1 = {max_colors}")
    print(f"So, the maximum number of colours is {max_colors}.")

solve_coloring_problem()