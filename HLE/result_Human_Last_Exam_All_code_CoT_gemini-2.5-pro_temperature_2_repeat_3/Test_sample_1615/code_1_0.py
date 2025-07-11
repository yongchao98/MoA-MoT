def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors for a non-complete graph with n vertices.
    """
    # The number of vertices in the graph G.
    n = 12345

    # Step 1: Explain the general theory.
    # The chromatic number of any graph with n vertices, chi(G), is at most n.
    # chi(G) = n if and only if G is a complete graph (K_n), where every vertex
    # is connected to every other vertex.
    # The problem states that G is NOT a complete graph.

    # Step 2: Establish the upper bound.
    # Since G is not a complete graph, there must be at least one pair of vertices
    # that are not adjacent. These two vertices can be assigned the same color.
    # This implies that a proper coloring with at most n-1 colors is always possible.
    # Therefore, the chromatic number must be less than n: chi(G) <= n - 1.

    # Step 3: Show the upper bound is achievable (tight).
    # To show that n-1 is the maximum, we need to find a non-complete graph on n
    # vertices that requires exactly n-1 colors.
    # Consider a graph formed by a complete subgraph on n-1 vertices (a K_{n-1})
    # and one additional isolated vertex.
    # This graph has (n-1) + 1 = n vertices. It is not complete because the
    # isolated vertex is not connected to the others.
    # The K_{n-1} part, by definition, requires n-1 distinct colors. So, the
    # chromatic number of the whole graph is at least n-1.

    # Step 4: Conclude the maximum number of colors.
    # From Step 2 and 3, we know chi(G) <= n-1, and there is a graph for which
    # chi(G) = n-1. So, the maximum is n-1.
    max_colors = n - 1

    # Step 5: Output the result and the equation.
    print(f"The number of vertices is n = {n}.")
    print("For any graph G that is not a complete graph, the chromatic number chi(G) is at most n - 1.")
    print("We can construct a non-complete graph with a chromatic number of n - 1 (e.g., K_{n-1} plus an isolated vertex).")
    print("Therefore, the maximum number of colors needed is n - 1.")
    print(f"The final calculation is: {n} - 1 = {max_colors}")

solve_graph_coloring_problem()