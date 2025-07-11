def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with n vertices, given that the graph is not complete.
    """
    
    # Number of vertices in the graph
    n_vertices = 12345
    
    # For any graph G with n vertices, the chromatic number chi(G) is n
    # if and only if G is the complete graph K_n.
    # The problem states that G is not a complete graph. This means there is at least
    # one pair of non-adjacent vertices. These two vertices can share the same color.
    # Therefore, the chromatic number must be less than n. The maximum possible
    # value is n - 1.
    
    # We can construct a graph that achieves this maximum. For example, a graph
    # consisting of a complete subgraph on (n-1) vertices (K_{n-1}) and one
    # additional isolated vertex. This graph has n vertices, is not complete,
    # and its chromatic number is n-1.
    
    # So, the maximum number of colors needed is n - 1.
    max_colors = n_vertices - 1
    
    print(f"The number of vertices in the graph G is n = {n_vertices}.")
    print("The graph G is not a complete graph.")
    print("The maximum number of colors needed for a proper vertex coloring of G is n - 1.")
    print("The final equation is:")
    print(f"{n_vertices} - 1 = {max_colors}")

solve_coloring_problem()