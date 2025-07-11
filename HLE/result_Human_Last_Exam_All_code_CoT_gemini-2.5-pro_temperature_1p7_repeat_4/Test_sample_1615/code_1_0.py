# The number of vertices in the graph G
n = 12345

# For any graph G with n vertices that is not a complete graph, the chromatic number chi(G)
# is at most n-1. This is because a non-complete graph must have at least one pair of
# non-adjacent vertices, which can be assigned the same color.
#
# To show that this maximum is achievable, consider a graph G' which is a complete graph K_n
# with one edge removed. This graph is not complete, but it contains a clique of size n-1.
# A clique of size n-1 requires n-1 colors.
#
# Therefore, the maximum number of colors needed for a proper coloring of a non-complete
# graph with n vertices is n-1.

# Calculate the result
max_colors = n - 1

# Print the final equation and the result
print(f"The number of vertices is n = {n}.")
print(f"The maximum number of colors required for a non-complete graph with n vertices is n - 1.")
print(f"Calculation: {n} - 1 = {max_colors}")