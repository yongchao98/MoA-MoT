import sys

def solve():
    """
    This function calculates the maximum number of colors needed to properly
    color a graph G with 12345 vertices, given that G is not a complete graph.
    """
    # The number of vertices in the graph G
    n = 12345

    # In graph theory, the chromatic number of a graph with n vertices is n
    # if and only if the graph is a complete graph (K_n).
    # A complete graph is one where every pair of distinct vertices is connected by an edge.

    # The problem specifies that our graph G is NOT the complete graph.
    # Therefore, its chromatic number must be strictly less than n.
    # The highest possible integer value less than n is n - 1.

    # We must confirm that a chromatic number of n-1 is achievable for a
    # non-complete graph with n vertices.
    # Consider a graph formed by a complete subgraph on (n-1) vertices (a K_{n-1})
    # and one additional isolated vertex. This graph has n vertices, is not complete,
    # and requires n-1 colors for a proper coloring.

    # Thus, the maximum number of colors that might be needed is n - 1.
    max_colors = n - 1

    print(f"The number of vertices is n = {n}.")
    print("For any graph with n vertices that is not a complete graph, the maximum number of colors required for a proper vertex coloring is n - 1.")
    print("The final calculation is:")
    print(f"{n} - 1 = {max_colors}")

solve()