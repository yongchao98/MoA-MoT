def solve():
    """
    Calculates the maximum number of colors needed for a non-complete graph with 12345 vertices.

    Let n be the number of vertices in a graph G.
    A key result in graph theory states that the chromatic number chi(G) is equal to n
    if and only if G is a complete graph (K_n).

    If G is not a complete graph, then its chromatic number chi(G) is at most n-1.
    The problem specifies the graph G has n=12345 vertices and is not a complete graph.
    Thus, the maximum possible chromatic number for G is n-1.

    We need to show this maximum is achievable. A graph consisting of a complete subgraph K_{n-1}
    and one isolated vertex is not a complete graph, and its chromatic number is n-1.
    So, the maximum is indeed n-1.
    """
    num_vertices = 12345

    # The maximum number of colors for a non-complete graph with n vertices is n - 1.
    max_colors = num_vertices - 1

    print(f"The number of vertices is n = {num_vertices}.")
    print("For a graph G with n vertices that is not a complete graph, the maximum number of colors needed is n - 1.")
    print(f"So, the maximum number of colors is {num_vertices} - 1 = {max_colors}.")

solve()