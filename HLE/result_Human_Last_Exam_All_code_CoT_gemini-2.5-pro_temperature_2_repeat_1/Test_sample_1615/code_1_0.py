def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph with 12345
    vertices, given that it is not a complete graph.
    """
    
    # The number of vertices in the graph G
    n = 12345

    # In graph theory, the chromatic number of a graph G, denoted Chi(G), is the minimum
    # number of colors needed for a proper vertex coloring.

    # A fundamental theorem states that for any graph with n vertices, Chi(G) <= n.
    # The equality, Chi(G) = n, holds if and only if the graph G is a complete graph (K_n),
    # where every vertex is connected to every other vertex.

    # The problem specifies that the graph G is not a complete graph.
    # Therefore, its chromatic number must be strictly less than n.
    # Chi(G) < n  which means Chi(G) <= n - 1.

    # This means the maximum possible number of colors is at most n - 1.
    
    # To confirm that n-1 is the maximum, we need to show that there exists a graph
    # with n vertices, which is not complete, but requires n-1 colors.
    # Consider a graph constructed by taking a complete graph on n-1 vertices (K_{n-1})
    # and adding one more vertex as an isolated point (not connected to anything).
    #
    # 1. Total vertices = (n - 1) + 1 = n.
    # 2. The graph is not a complete K_n because the isolated vertex is not connected to others.
    # 3. The K_{n-1} part of the graph requires n-1 colors to be colored properly.
    # 4. The single isolated vertex can be colored with any of the existing n-1 colors.
    #
    # The chromatic number of this constructed graph is therefore n-1.

    # Since the maximum number of colors is at most n-1, and we have found a graph
    # that requires exactly n-1 colors, the maximum is n-1.

    max_colors = n - 1

    print(f"Let n be the number of vertices in the graph G. Here, n = {n}.")
    print("The problem states that G is not a complete graph.")
    print("A graph on n vertices requires n colors if and only if it is the complete graph K_n.")
    print("Since G is not a complete graph, the number of colors must be less than n.")
    print(f"This means the maximum number of colors is at most n - 1, which is {n} - 1 = {max_colors}.")
    print("\nWe can construct a graph with a chromatic number of n - 1.")
    print(f"Consider a graph with a complete subgraph on {n-1} vertices and one extra, isolated vertex.")
    print(f"This graph has {n} vertices, is not complete, and requires {n-1} colors.")
    print("\nTherefore, the maximum number of colors needed is given by the equation:")
    print(f"{n} - 1 = {max_colors}")

solve_graph_coloring_problem()