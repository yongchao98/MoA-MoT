def solve_graph_coloring_problem():
    """
    This function calculates the maximum number of colors needed for a proper vertex
    coloring of a graph with 12345 vertices that is not a complete graph.

    The logic is as follows:
    1. A graph with n vertices has a chromatic number of n if and only if it is a
       complete graph (K_n).
    2. Since the graph G is not a complete graph, its chromatic number must be
       less than n, i.e., chi(G) <= n - 1.
    3. We can construct a graph with n vertices that is not complete and has a
       chromatic number of n-1. An example is a complete graph on n-1 vertices
       plus one isolated vertex.
    4. Therefore, the maximum possible chromatic number is n - 1.
    """
    
    num_vertices = 12345
    
    # The maximum number of colors is one less than the number of vertices.
    max_colors = num_vertices - 1
    
    # Print the equation as requested.
    print(f"The number of vertices is n = {num_vertices}.")
    print("The maximum number of colors needed for a non-complete graph with n vertices is n - 1.")
    print(f"Calculation: {num_vertices} - 1 = {max_colors}")

solve_graph_coloring_problem()