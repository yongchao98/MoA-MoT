# The number of vertices in the graph G.
n = 12345

# The problem asks for the maximum number of colors required for a proper vertex coloring
# of a graph G with n vertices, given that G is not a complete graph.

# 1. A complete graph with n vertices, K_n, is the only graph that requires n colors.
#    In a complete graph, every vertex is connected to every other vertex.

# 2. Since the graph G is not a complete graph, it means there is at least one pair of
#    vertices that are not connected by an edge. This implies its chromatic number,
#    chi(G), must be less than n. Thus, chi(G) <= n - 1.

# 3. To find the maximum number of colors, we need to check if it's possible for a
#    non-complete graph to require n - 1 colors. Consider a graph G' that is a
#    complete graph K_n with just one edge removed. This graph is not complete,
#    and it can be shown that its chromatic number is exactly n - 1.

# 4. Therefore, the maximum number of colors needed for a non-complete graph with
#    n vertices is n - 1.

# Calculate the result for n = 12345.
max_colors = n - 1

# Print the final equation showing the calculation.
print(f"For a non-complete graph with {n} vertices, the maximum number of colors required for a proper coloring is {n} - 1.")
print(f"The final equation is: {n} - 1 = {max_colors}")