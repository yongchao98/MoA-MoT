# The number of vertices in the graph G.
n = 12345

# The problem asks for the maximum number of colors needed for a proper vertex
# coloring of a graph with n vertices, given that the graph is not a complete graph.

# Based on graph theory, the chromatic number of any graph with n vertices
# is at most n. The equality holds if and only if the graph is a complete graph (K_n).
# Since our graph G is not a complete graph, its chromatic number must be at most n - 1.

# To show that n - 1 is the highest possible value, we can construct a graph that
# requires n - 1 colors. A simple example is a graph consisting of a complete
# subgraph on n-1 vertices (K_{n-1}) and one additional isolated vertex. This graph
# has n vertices, is not complete, and requires n-1 colors.

# Thus, the maximum number of colors required is n - 1.
max_colors = n - 1

print(f"Let n be the number of vertices, where n = {n}.")
print("The maximum number of colors needed to properly color a non-complete graph with n vertices is given by the equation: n - 1.")
print("The result is:")
print(f"{n} - 1 = {max_colors}")