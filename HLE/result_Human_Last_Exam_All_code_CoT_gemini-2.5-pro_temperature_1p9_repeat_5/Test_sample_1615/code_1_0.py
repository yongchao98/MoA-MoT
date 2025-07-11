# The number of vertices in the graph
n = 12345

# According to graph theory, the maximum number of colors needed to properly
# color a graph G with n vertices is n if and only if G is the complete graph K_n.
# If the graph G is not a complete graph, its chromatic number is at most n-1.
# This upper bound is achieved by a graph constructed by removing a single edge
# from the complete graph K_n. Such a graph contains a clique of size n-1,
# and therefore requires at least n-1 colors.
# Thus, the maximum number of colors is exactly n-1.

# Calculate the maximum number of colors
max_colors = n - 1

# Print the final equation as requested.
print(f"The number of vertices is n = {n}.")
print(f"Since the graph is not complete, the maximum number of colors required is n - 1.")
print(f"The final calculation is: {n} - 1 = {max_colors}")