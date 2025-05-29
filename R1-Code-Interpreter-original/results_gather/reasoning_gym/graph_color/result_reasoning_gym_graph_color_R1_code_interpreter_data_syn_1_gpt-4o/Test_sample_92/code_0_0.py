# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (1, 2)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Assign colors to the connected vertices
for edge in edges:
    v1, v2 = edge
    # Assign the first available color to v1
    if v1 not in color_map:
        color_map[v1] = colors[0]
    # Assign a different color to v2
    for color in colors:
        if color != color_map[v1]:
            color_map[v2] = color
            break

# Assign colors to the remaining vertices
for vertex in vertices:
    if vertex not in color_map:
        # Assign the first available color
        color_map[vertex] = colors[0]

# Print the color map
print(color_map)