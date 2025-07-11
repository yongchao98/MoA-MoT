def solve_graph_coloring_problem():
    """
    This function calculates the maximum number of colors needed to properly
    color a graph G with n vertices, given that G is not a complete graph.
    """
    
    # The number of vertices in the graph G.
    num_vertices = 12345
    
    # According to graph theory, a graph with n vertices has a chromatic number of n
    # if and only if it is the complete graph K_n.
    # Since the graph G is not a complete graph, its chromatic number must be less than n.
    # So, the maximum possible chromatic number is n - 1.
    # This maximum is achieved, for example, by a graph composed of a complete
    # subgraph K_{n-1} and one isolated vertex.
    
    max_colors = num_vertices - 1
    
    # Print the equation as requested.
    print(f"{num_vertices} - 1 = {max_colors}")

solve_graph_coloring_problem()