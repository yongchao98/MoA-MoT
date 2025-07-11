def solve():
    """
    Calculates the maximum number of colours needed to properly colour the vertices
    of a graph with 12345 vertices that is not a complete graph.
    """
    num_vertices = 12345
    
    # For a graph G with n vertices, chi(G) <= n.
    # Equality holds if and only if G is a complete graph (K_n).
    # Since the graph is not complete, its chromatic number must be less than n.
    # So, chi(G) <= n - 1.
    # We can construct a graph that achieves this bound. For example, a graph
    # consisting of a complete graph on n-1 vertices and one isolated vertex
    # is not a complete graph and has a chromatic number of n-1.
    # Therefore, the maximum number of colors is n - 1.
    
    max_colors = num_vertices - 1
    
    # Print the equation as requested.
    print(f"{num_vertices} - 1 = {max_colors}")

solve()