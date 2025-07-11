def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper coloring of a
    graph with n vertices, given that the graph is not complete.

    The chromatic number chi(G) of a graph G with n vertices is n if and only if
    G is a complete graph (K_n).
    Since the graph G is not a complete graph, its chromatic number must be
    less than n. The maximum possible value is n-1.

    This can be achieved by a graph consisting of a complete subgraph on n-1
    vertices (a K_{n-1} clique) and one additional isolated vertex. Such a
    graph has n vertices, is not complete, and its chromatic number is n-1.
    """
    num_vertices = 12345
    max_colors = num_vertices - 1

    print(f"Let n be the number of vertices, n = {num_vertices}.")
    print("The graph is not a complete graph.")
    print("The maximum number of colors required is n - 1.")
    print(f"The calculation is: {num_vertices} - 1 = {max_colors}")

solve_coloring_problem()