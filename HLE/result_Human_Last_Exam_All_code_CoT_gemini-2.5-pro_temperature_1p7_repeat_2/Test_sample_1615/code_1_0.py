def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper coloring of a
    graph with n vertices that is not a complete graph.
    """
    
    # Number of vertices in the graph G
    n = 12345
    
    # According to graph theory, a graph G with n vertices has a chromatic number
    # chi(G) = n if and only if G is a complete graph (K_n).
    # Since the problem states G is not a complete graph, its chromatic number
    # must be less than n. The maximum possible value is n - 1.
    # This value is achieved by a complete graph with one edge removed (K_n - e),
    # which contains a clique of size n-1.
    
    max_colors = n - 1
    
    # Print the equation as requested
    print(f"The maximum number of colors is determined by the equation: {n} - 1 = {max_colors}")

solve_coloring_problem()