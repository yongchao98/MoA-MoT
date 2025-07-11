def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring of
    a graph with 12345 vertices that is not a complete graph.
    """
    # The number of vertices in the graph G
    n = 12345

    # According to Brooks' Theorem and properties of graph coloring:
    # 1. For any graph G with n vertices, the chromatic number chi(G) <= n.
    # 2. The equality chi(G) = n holds if and only if G is a complete graph (K_n).
    # 3. Since the graph G is not a complete graph, its chromatic number must be less than n.
    #    Therefore, chi(G) <= n - 1.

    # To show that n - 1 is the maximum possible value, we can construct a graph H
    # with n vertices that is not complete and has a chromatic number of n - 1.
    # Let H be the graph formed by a complete subgraph on n - 1 vertices (K_{n-1})
    # and one additional isolated vertex.
    # This graph H has n vertices, is not a complete graph, and its chromatic number
    # is chi(K_{n-1}) = n - 1.

    # Thus, the maximum number of colors needed is n - 1.
    max_colors = n - 1

    print(f"The graph G has n = {n} vertices and is not a complete graph.")
    print("The maximum number of colors required for a proper coloring of such a graph is n - 1.")
    print(f"Therefore, the maximum number of colours is {n} - 1 = {max_colors}.")

solve_graph_coloring_problem()