def solve_coloring_problem():
    """
    Calculates the maximum number of colors for a non-complete graph with n vertices.
    """
    # The number of vertices in the graph G
    n = 12345
    
    # The chromatic number of a graph with n vertices is at most n.
    # Equality holds if and only if the graph is a complete graph (K_n).
    # Since the graph G is not a complete graph, its chromatic number must be less than n.
    # So, the maximum possible chromatic number is n - 1.
    # This can be achieved by a graph consisting of a K_{n-1} clique and one isolated vertex.
    max_colors = n - 1
    
    print(f"Let n be the number of vertices in the graph, where n = {n}.")
    print("The graph is not a complete graph.")
    print("The maximum number of colors needed for a proper vertex coloring of such a graph is n - 1.")
    print("The calculation is:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()