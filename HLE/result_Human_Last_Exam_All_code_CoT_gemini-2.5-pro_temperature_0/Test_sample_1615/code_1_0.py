def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph
    with 12345 vertices that is not a complete graph.
    """
    
    # Number of vertices in the graph G
    n = 12345
    
    # For any graph with n vertices, the chromatic number chi(G) is at most n.
    # The equality chi(G) = n holds if and only if G is a complete graph (K_n).
    # The problem states that G is not a complete graph.
    # Therefore, the chromatic number of G must be strictly less than n.
    # This means chi(G) <= n - 1.
    
    # To show this is the maximum, we must show it's achievable.
    # Consider a graph formed by a complete subgraph on (n-1) vertices (a K_{n-1})
    # and one additional isolated vertex. This graph has n vertices, is not complete,
    # and its chromatic number is n-1.
    
    # Therefore, the maximum number of colors is n - 1.
    max_colors = n - 1
    
    print(f"Let n be the number of vertices in the graph G. Here, n = {n}.")
    print("The problem states that G is not a complete graph.")
    print("The chromatic number of any graph with n vertices, chi(G), is at most n.")
    print("The case chi(G) = n occurs only when G is a complete graph, K_n.")
    print("Since G is not a complete graph, its chromatic number must be less than n.")
    print(f"This gives an upper bound: chi(G) <= n - 1.")
    print(f"The maximum number of colors is therefore at most {n} - 1 = {max_colors}.")
    print(f"This maximum is achievable for a graph like K_{n-1} with an isolated vertex.")
    print(f"Thus, the maximum number of colours needed is {max_colors}.")

solve_coloring_problem()