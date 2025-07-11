def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph G with n vertices, given that G is not a complete graph.
    """
    # Number of vertices in the graph G
    n = 12345

    # According to graph theory, the chromatic number of a graph with n vertices, chi(G),
    # is equal to n if and only if G is a complete graph (K_n).
    # Since the problem states G is not a complete graph, its chromatic number must be less than n.
    # Therefore, the maximum possible chromatic number is n - 1.
    # We can prove this is achievable by constructing a graph, for example, a K_{n-1}
    # subgraph with one additional isolated vertex. This graph has n vertices, is not complete,
    # and has a chromatic number of n - 1.

    # Perform the calculation
    max_colors = n - 1
    
    # Output the explanation and the final equation.
    print(f"The number of vertices in the graph is n = {n}.")
    print("For any graph G with n vertices that is not a complete graph, the chromatic number chi(G) is at most n-1.")
    print(f"The maximum number of colors needed is therefore n - 1.")
    print("The final equation is:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()
