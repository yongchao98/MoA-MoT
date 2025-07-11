def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph
    with 12345 vertices that is not a complete graph.
    """
    
    # The number of vertices in the graph G.
    num_vertices = 12345
    
    # For any graph G with n vertices, the chromatic number chi(G) is at most n.
    # The equality chi(G) = n holds if and only if G is the complete graph K_n.
    # The problem states that G is not a complete graph.
    # Therefore, the chromatic number of G must be strictly less than n.
    # The maximum possible value for chi(G) is n - 1.
    
    # This maximum is achievable. For example, a graph consisting of a
    # complete subgraph K_{n-1} and one additional isolated vertex has n vertices,
    # is not complete, and its chromatic number is n-1.
    
    # So, the maximum number of colors needed is n - 1.
    max_colors = num_vertices - 1
    
    # Print the equation showing the calculation.
    print(f"The maximum number of colours is given by the equation:")
    print(f"{num_vertices} - 1 = {max_colors}")

solve_graph_coloring_problem()