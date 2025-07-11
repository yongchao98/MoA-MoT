def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with n vertices, given it is not a complete graph.
    """
    
    # Number of vertices in the graph G
    n = 12345
    
    # According to graph theory, the chromatic number chi(G) of any graph with n
    # vertices is at most n. The equality chi(G) = n holds if and only if G is a
    # complete graph, K_n.
    
    # The problem states that G is not a complete graph.
    # Therefore, chi(G) must be strictly less than n.
    # So, chi(G) <= n - 1.
    
    # To show that n - 1 is the maximum possible number of colors, we need to show
    # that a graph exists with n vertices (not complete) that requires exactly n-1 colors.
    # Consider a graph formed by a complete graph on n-1 vertices (K_{n-1}) and one
    # isolated vertex. This graph has n vertices, is not complete, and its
    # chromatic number is n-1.
    
    # Thus, the maximum number of colors needed is n - 1.
    max_colors = n - 1
    
    print(f"Let n be the number of vertices, so n = {n}.")
    print("The graph is not a complete graph, so its chromatic number chi(G) < n.")
    print("This means the maximum number of colors is at most n - 1.")
    print("An example requiring n - 1 colors is a complete graph K_{n-1} plus one isolated vertex.")
    print("Therefore, the maximum number of colors required is given by the equation:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()
