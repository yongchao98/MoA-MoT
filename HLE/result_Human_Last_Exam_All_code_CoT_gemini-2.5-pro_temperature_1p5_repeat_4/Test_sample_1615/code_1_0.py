def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with 12345 vertices that is not a complete graph.
    """
    
    # The number of vertices in the graph G.
    n = 12345
    
    # A graph on n vertices requires n colors if and only if it is the complete graph Kn.
    # In a complete graph, every vertex is connected to every other vertex, so each
    # vertex must have a unique color.
    
    # The problem states that G is NOT a complete graph.
    # Therefore, its chromatic number must be less than n.
    # The maximum possible chromatic number is thus n - 1.
    
    # We can prove this is achievable by constructing a graph G' that is not complete
    # and has a chromatic number of n - 1. For instance, consider a graph formed by
    # a complete subgraph on n-1 vertices (K_{n-1}) and one isolated vertex.
    # This graph has n vertices, is not complete, but requires n-1 colors for the K_{n-1} part.
    # Thus, the maximum number of colors required is n - 1.
    
    max_colors = n - 1
    
    # Print the equation as requested.
    print(f"The number of vertices is n = {n}.")
    print("For a graph G with n vertices, the chromatic number is n if and only if G is a complete graph.")
    print("Since G is not a complete graph, its chromatic number must be at most n-1.")
    print("This maximum is achievable, for example, by a graph consisting of a K_{n-1} clique and one isolated vertex.")
    print("Therefore, the maximum number of colours is:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()