# The number of vertices in the graph G
num_vertices = 12345

# According to graph theory, the maximum chromatic number for a graph with n vertices
# that is not a complete graph is n-1.
# A complete graph Kn requires n colors. Any graph other than Kn is missing at least
# one edge. Let u and v be two vertices that are not connected. We can assign the
# same color to u and v. All other n-2 vertices can be assigned their own unique
# colors. This gives a valid coloring with at most 1 + (n-2) = n-1 colors.
# This maximum is achievable, for example, by a graph consisting of a complete graph
# on n-1 vertices (K_{n-1}) and one isolated vertex. Such a graph requires n-1 colors.

max_colors = num_vertices - 1

# Print the equation showing the calculation.
print(f"The maximum number of colours required is given by the equation:")
print(f"{num_vertices} - 1 = {max_colors}")