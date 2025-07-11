# The number of vertices in the graph G
num_vertices = 12345

# For any graph G with n vertices, the chromatic number chi(G) is at most n.
# The chromatic number is n if and only if the graph is a complete graph (K_n).
# Since the graph G is not a complete graph, its chromatic number must be less than n.
# So, the maximum possible chromatic number is n - 1.
# We can construct a graph with n-1 vertices forming a clique (K_{n-1}) and one isolated vertex.
# This graph has n vertices, is not complete, and has a chromatic number of n - 1.
max_colors = num_vertices - 1

# Print the final calculation as an equation.
print(f"{num_vertices} - 1 = {max_colors}")