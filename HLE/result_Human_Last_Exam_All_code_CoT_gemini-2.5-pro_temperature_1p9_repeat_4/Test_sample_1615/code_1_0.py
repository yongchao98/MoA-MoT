def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph
    with 12345 vertices that is not a complete graph.
    """
    num_vertices = 12345

    # A graph G with n vertices has a chromatic number of n if and only if
    # G is a complete graph (K_n).
    # Since the graph in the problem is not a complete graph, its chromatic
    # number must be less than n.
    max_colors = num_vertices - 1

    print("The problem asks for the maximum number of colors needed to color a graph G with 12345 vertices, given it is not a complete graph.")
    print("Let n be the number of vertices.")
    print(f"Here, n = {num_vertices}.")
    print("\nA graph on n vertices has a chromatic number of n if and only if it is a complete graph.")
    print("Since G is not a complete graph, its chromatic number must be strictly less than n.")
    print(f"This implies the number of colors must be at most n - 1.")
    print("\nTo confirm this is the maximum, we need to show that a graph needing n - 1 colors exists.")
    print("Consider a graph formed by a complete subgraph on n - 1 vertices (a (n-1)-clique) and one other vertex not connected to all vertices of the clique.")
    print("Such a graph is not complete, but the (n-1)-clique itself requires n - 1 distinct colors.")
    print("Therefore, the maximum number of colors required is indeed n - 1.")
    print("\nThe final calculation is:")
    print(f"{num_vertices} - 1 = {max_colors}")

solve_coloring_problem()