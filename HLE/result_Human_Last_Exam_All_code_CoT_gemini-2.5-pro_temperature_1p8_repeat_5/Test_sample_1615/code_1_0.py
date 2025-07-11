def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph with n vertices,
    given that the graph is not a complete graph.
    """
    # The number of vertices in the graph G.
    n = 12345

    # The problem is to find the maximum chromatic number chi(G) for a graph G with
    # n vertices, given that G is not the complete graph K_n.

    # 1. Finding an upper bound:
    # A standard result in graph theory is that for any graph G with n vertices,
    # its chromatic number chi(G) is less than or equal to n.
    # The equality chi(G) = n holds if and only if G is the complete graph K_n.
    # Since we are given that G is *not* a complete graph, its chromatic number must
    # be strictly less than n. Therefore, chi(G) <= n - 1.

    # 2. Showing the upper bound is achievable:
    # To show that n-1 is the maximum possible value, we need to find an example
    # of a graph that is not complete and requires n-1 colors.
    # Consider a graph G* constructed from a complete graph on n-1 vertices (K_{n-1})
    # and a single additional vertex that is isolated.
    # - The total number of vertices in G* is (n-1) + 1 = n.
    # - G* is not a complete graph because the isolated vertex is not connected to others.
    # - To color the K_{n-1} part of G*, we need n-1 distinct colors. The isolated
    #   vertex can take any of these colors. So, the chromatic number of G* is n-1.

    # Conclusion:
    # Since chi(G) <= n-1 for any such graph and there exists a graph that achieves
    # this bound, the maximum number of colors required is n - 1.
    max_colors = n - 1

    # Print the final equation as requested.
    print(f"The number of vertices is n = {n}.")
    print("The maximum number of colours for a proper vertex colouring of a graph G with n vertices, where G is not the complete graph, is n - 1.")
    print(f"Calculation: {n} - 1 = {max_colors}")

solve_coloring_problem()