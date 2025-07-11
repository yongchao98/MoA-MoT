def solve_max_colors():
    """
    This function determines the maximum number of colors needed for a proper vertex coloring
    of a graph G with 12345 vertices, given that G is not a complete graph.
    """
    
    # Number of vertices in the graph
    num_vertices = 12345

    # According to graph theory, a graph with n vertices requires n colors
    # for a proper vertex coloring if and only if it is the complete graph K_n.
    # In a complete graph, every vertex is connected to every other vertex.

    # The problem states that the graph G is NOT a complete graph.
    # This implies that its chromatic number, chi(G), must be less than the number of vertices.
    # So, chi(G) < 12345.
    
    # The maximum integer value less than 12345 is 12344.
    # This means the maximum possible number of colors is at most 12344.
    max_colors_upper_bound = num_vertices - 1

    # We must confirm if this maximum is actually achievable.
    # We can construct a graph with 12345 vertices that is not complete and
    # requires 12344 colors. For example, consider a graph made of two components:
    # 1. A complete graph K_12344 (12344 vertices, all connected to each other).
    # 2. A single isolated vertex.
    # This graph has 12345 vertices in total and is not complete (due to the isolated vertex).
    # The K_12344 component requires 12344 colors. The isolated vertex can reuse any of these colors.
    # Thus, the chromatic number of this graph is 12344.
    
    # Since the chromatic number cannot be greater than 12344, and we have found a
    # graph for which it is exactly 12344, this is the maximum possible value.
    
    final_answer = max_colors_upper_bound

    print("The number of vertices is n = {}.".format(num_vertices))
    print("A graph requires n colors if and only if it is a complete graph.")
    print("Since the graph is not complete, the number of colors must be less than n.")
    print("The maximum number of colors is therefore n - 1.")
    print("The equation is: {} - 1 = {}".format(num_vertices, final_answer))

solve_max_colors()