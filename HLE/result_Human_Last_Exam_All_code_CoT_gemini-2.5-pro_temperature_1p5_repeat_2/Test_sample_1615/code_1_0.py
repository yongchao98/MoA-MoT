# The number of vertices in the graph G
n = 12345

# For a graph G with n vertices, the chromatic number is n if and only if G is a complete graph.
# Since the graph is not complete, its chromatic number must be less than n.
# The maximum possible chromatic number is n-1.
# This can be achieved by a graph consisting of a complete subgraph K_{n-1} and one isolated vertex.
max_colors = n - 1

# Print the equation and the final answer.
print(f"{n} - 1 = {max_colors}")