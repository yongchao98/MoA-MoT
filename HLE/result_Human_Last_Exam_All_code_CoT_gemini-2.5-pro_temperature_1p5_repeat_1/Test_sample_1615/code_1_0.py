def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph G with n vertices, given that G is not a complete graph.
    """
    
    # The number of vertices in the graph G.
    n = 12345

    # Step 1: Establish an upper bound for the chromatic number, chi(G).
    # A graph with n vertices requires at most n colors. This happens when the graph
    # is a complete graph (K_n), where every vertex needs a unique color.
    # The problem states that G is NOT a complete graph. This implies there is at
    # least one pair of vertices, say u and v, that are not connected by an edge.
    # Since u and v are not adjacent, they can legally be assigned the same color.
    #
    # A valid coloring with n-1 colors can always be constructed. For example:
    # Color the n-2 vertices in V(G) \ {u, v} with n-2 unique colors.
    # Then, assign a new, (n-1)-th color to both u and v.
    # This is a proper coloring using n-1 colors. Therefore, the minimum number of
    # colors needed, chi(G), must be less than or equal to n-1.
    # So, chi(G) <= n - 1.

    # Step 2: Show that this upper bound is achievable.
    # We need to construct a graph G on n vertices which is not complete and has
    # a chromatic number of exactly n-1.
    # Consider the graph G formed by the disjoint union of a complete graph on n-1
    # vertices (K_{n-1}) and a single isolated vertex (K_1).
    # - This graph has (n-1) + 1 = n vertices.
    # - It is not a complete graph because the isolated vertex is not connected to any
    #   other vertices.
    # - The subgraph K_{n-1} requires n-1 distinct colors for a proper coloring.
    #   Thus, the chromatic number of the entire graph G must be at least n-1.
    #   So, chi(G) >= n - 1.

    # Step 3: Conclude and calculate.
    # From Step 1 (chi(G) <= n-1) and Step 2 (chi(G) >= n-1 for some G), we can
    # conclude that the maximum number of colors required is exactly n-1.
    
    maximum_colors = n - 1
    
    # Print the final result as an equation.
    print(f"The number of vertices is n = {n}.")
    print("For a graph G on n vertices that is not complete, the chromatic number chi(G) is at most n-1.")
    print("This maximum is achieved by a graph such as a K_{n-1} clique with one isolated vertex.")
    print("Therefore, the maximum number of colours needed is n - 1.")
    print(f"Calculation: {n} - 1 = {maximum_colors}")

solve_coloring_problem()