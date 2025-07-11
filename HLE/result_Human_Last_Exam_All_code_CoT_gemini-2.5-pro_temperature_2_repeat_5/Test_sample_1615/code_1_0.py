def solve_coloring_problem():
    """
    Calculates the maximum number of colors for a non-complete graph with a given number of vertices.
    """
    num_vertices = 12345
    
    # A graph G with n vertices has a chromatic number of n if and only if it is a complete graph (K_n).
    # The problem states that the graph G is not a complete graph.
    # Therefore, its chromatic number must be less than the number of vertices.
    # chi(G) <= n - 1
    
    max_colors = num_vertices - 1
    
    # To show that this maximum is achievable, consider a graph G constructed as follows:
    # Let G be the disjoint union of a complete graph on (n-1) vertices (K_{n-1}) and a single isolated vertex.
    # This graph has n vertices.
    # It is not a complete graph.
    # Its chromatic number is n-1, as the K_{n-1} component requires n-1 colors.
    
    print(f"The number of vertices in the graph is n = {num_vertices}.")
    print("The graph is not a complete graph.")
    print("The maximum number of colors needed for a proper vertex coloring of such a graph is n - 1.")
    print("This is because any graph G on n vertices has a chromatic number of n if and only if G is a complete graph.")
    print("Since G is not complete, its chromatic number must be less than n.")
    print(f"The maximum possible chromatic number is {num_vertices - 1}.")
    print("\nThe final calculation is:")
    print(f"{num_vertices} - 1 = {max_colors}")

solve_coloring_problem()
