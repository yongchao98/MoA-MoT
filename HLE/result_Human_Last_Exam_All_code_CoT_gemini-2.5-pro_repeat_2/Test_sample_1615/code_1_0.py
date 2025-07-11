# The number of vertices in the graph G
num_vertices = 12345

# The problem asks for the maximum number of colors needed to properly color a graph
# with 'num_vertices' vertices, given that the graph is not a complete graph.

# The chromatic number of a complete graph Kn is n.
# If a graph G is not complete, it means there's at least one pair of vertices
# that are not adjacent. These two vertices can be given the same color.
# This implies that the chromatic number of a non-complete graph on n vertices
# can be at most n-1.

# To show this maximum is achievable, consider a graph formed by a complete subgraph
# on (n-1) vertices (a K_{n-1}) and one isolated vertex. This graph is not complete,
# has n vertices, and its chromatic number is n-1.

# Therefore, the maximum number of colors is num_vertices - 1.
max_colors = num_vertices - 1

# The final code should output each number in the final equation.
print(f"The maximum number of colors for a non-complete graph with {num_vertices} vertices is:")
print(f"{num_vertices} - 1 = {max_colors}")