# Initialize the graph coloring
vertices = list(range(10))
edges = [(6, 8), (8, 9)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Color the vertices connected by edges
color_map[6] = 1
color_map[8] = 2
color_map[9] = 1

# Color the remaining vertices
for vertex in vertices:
    if vertex not in color_map:
        color_map[vertex] = 1

# Print the color map
print(color_map)