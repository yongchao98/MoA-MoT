# Define the number of vertices in the graph G
num_vertices = 12345

# The chromatic number of a graph is the minimum number of colors needed for a proper vertex coloring.
# A graph with 'n' vertices has a chromatic number of at most 'n'.
# This maximum is only reached if the graph is a complete graph (K_n), where every vertex is connected to every other vertex.

# The problem states that the graph G is NOT a complete graph.
# If a graph is not complete, it means there exists at least one pair of vertices that are not connected by an edge.
# Let these two non-adjacent vertices be 'u' and 'v'.
# Because they are not connected, 'u' and 'v' can be assigned the same color in a proper coloring.

# This implies that we can always color a non-complete graph with at most n-1 colors.
# For example, by assigning one color to the non-adjacent pair {u, v} and a unique color to each of the remaining (n-2) vertices.
# This gives a valid coloring with 1 + (n-2) = n-1 colors.
# So, the maximum possible number of colors is at most n-1.

# To show that this maximum is achievable, consider a graph constructed by taking a complete graph on (n-1) vertices
# and adding one more vertex that is isolated (not connected to any other vertex).
# This graph has n vertices, is not complete, and requires n-1 colors due to the complete subgraph.
# Therefore, the maximum number of colors for a non-complete graph on n vertices is exactly n-1.

# Calculate the result for n = 12345
max_colors = num_vertices - 1

# Print the final equation as requested
print("The number of vertices is n = 12345.")
print("The graph is not a complete graph, so the maximum number of colors required is n - 1.")
print("The final equation is:")
print(f"{num_vertices} - 1 = {max_colors}")