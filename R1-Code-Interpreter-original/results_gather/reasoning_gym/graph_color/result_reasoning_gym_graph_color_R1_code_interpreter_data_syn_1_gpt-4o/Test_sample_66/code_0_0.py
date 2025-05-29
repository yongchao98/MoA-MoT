# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 2), (1, 6)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Assign colors to vertices connected by edges
color_map[0] = colors[0]
color_map[2] = colors[1]
color_map[1] = colors[0]
color_map[6] = colors[1]

# Assign colors to the remaining vertices
for vertex in vertices:
    if vertex not in color_map:
        color_map[vertex] = colors[0]  # Assign the first color available

# Print the color map
print(color_map)