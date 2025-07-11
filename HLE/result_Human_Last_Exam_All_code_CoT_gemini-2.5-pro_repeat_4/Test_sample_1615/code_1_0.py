def solve_coloring_problem():
    """
    Calculates the maximum number of colours needed to properly colour a graph
    with 12345 vertices, given that it is not a complete graph.
    """
    
    # The number of vertices in the graph G.
    num_vertices = 12345

    # The problem asks for the maximum number of colours required for a proper vertex colouring
    # of a graph G with 'num_vertices', given that G is not a complete graph.
    # This value is the maximum possible chromatic number, chi(G), for any such graph.

    # A key theorem in graph theory provides an upper bound for the chromatic number of a non-complete graph.
    # Theorem: If a graph G on n vertices is not a complete graph, then chi(G) <= n - 1.
    
    # Proof sketch:
    # Since G is not a complete graph, there must be at least one vertex, v, that is not
    # connected to all other n-1 vertices. Thus, the degree of v is at most n-2.
    # We can use a greedy colouring algorithm and order the vertices such that v is coloured last.
    # When colouring v, its neighbours (at most n-2 of them) have already been coloured.
    # These neighbours use at most n-2 distinct colours. So, from a palette of n-1 colours,
    # there is guaranteed to be at least one colour available for v.
    # This proves that any non-complete graph on n vertices can be coloured with at most n-1 colours.

    # So, for n = 12345, the maximum number of colours is at most n - 1.
    max_colors_upper_bound = num_vertices - 1

    # To show that this maximum is achievable, we need to find an example of a graph that
    # satisfies the conditions and requires exactly n - 1 colours.
    # Consider a graph G constructed as the union of two disjoint components:
    # 1. A complete graph on n-1 vertices (K_{n-1}).
    # 2. A single isolated vertex.
    
    # This graph G has (n-1) + 1 = n vertices.
    # It is not a complete graph because the isolated vertex is not connected to any other vertices.
    # The chromatic number of G is the maximum of the chromatic numbers of its components.
    # The K_{n-1} component requires exactly n-1 colours. The isolated vertex requires 1 colour.
    # Thus, chi(G) = max(n-1, 1) = n-1.

    # Since the maximum number of colours is at most n-1, and we have found a graph that requires
    # exactly n-1 colours, the maximum possible number of colours is n-1.

    max_colours_needed = num_vertices - 1

    print("The problem is to find the maximum number of colours to properly colour a graph G with n vertices, where n=12345 and G is not a complete graph.")
    print("Based on graph theory, the maximum number of colours required is n - 1.")
    print("\nThe final calculation is:")
    print(f"{num_vertices} - 1 = {max_colours_needed}")

solve_coloring_problem()