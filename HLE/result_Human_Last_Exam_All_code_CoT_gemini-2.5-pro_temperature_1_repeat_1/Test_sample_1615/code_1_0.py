def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with n vertices that is not a complete graph.
    """
    # The number of vertices in the graph G
    n = 12345

    # A fundamental theorem in graph theory states that the chromatic number chi(G) of a graph G
    # with n vertices is n if and only if G is a complete graph (K_n).
    # If G is not a complete graph, its chromatic number must be less than n.
    # So, the maximum possible chromatic number is n - 1.
    
    # We can construct a graph that achieves this maximum: a complete graph on n-1 vertices (K_{n-1})
    # plus one isolated vertex. This graph has n vertices, is not complete, and its chromatic number is n-1.

    # Therefore, the maximum number of colors required is n - 1.
    max_colors = n - 1

    print(f"The number of vertices in the graph is n = {n}.")
    print("Since the graph is not a complete graph, the maximum number of colors needed for a proper coloring is n - 1.")
    print("The final calculation is:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()