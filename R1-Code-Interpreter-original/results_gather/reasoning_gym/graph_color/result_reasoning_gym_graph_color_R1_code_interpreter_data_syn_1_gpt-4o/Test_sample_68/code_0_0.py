# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(5, 6), (7, 9)]
colors = [1, 2, 3]

# Initialize the coloring map
coloring = {}

# Assign colors to the vertices connected by edges
for edge in edges:
    v1, v2 = edge
    # Assign the first color to the first vertex if not already colored
    if v1 not in coloring:
        coloring[v1] = colors[0]
    # Assign a different color to the second vertex
    for color in colors:
        if color != coloring[v1]:
            coloring[v2] = color
            break

# Assign colors to the remaining vertices
for vertex in vertices:
    if vertex not in coloring:
        coloring[vertex] = colors[0]  # Assign the first color

# Print the resulting coloring
print(coloring)