# Number of vertices in the graph
num_vertices = 12345

# For any graph G with n vertices, the chromatic number chi(G) is at most n.
# chi(G) = n if and only if G is a complete graph.
# The problem states that the graph is NOT a complete graph.
# Therefore, the chromatic number must be less than n.
# The maximum possible chromatic number is n - 1.
# We can construct such a graph, for example, by taking a complete graph on n-1 vertices
# and adding one isolated vertex. This graph has n vertices, is not complete,
# and its chromatic number is n-1.

max_colors = num_vertices - 1

# Output the equation and the result
print(f"The maximum number of colors needed is given by the equation:")
print(f"{num_vertices} - 1 = {max_colors}")