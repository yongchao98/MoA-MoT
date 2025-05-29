# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 5)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Assign colors
color_map[1] = 1
color_map[5] = 2

# Assign color 1 to all other vertices
for vertex in vertices:
    if vertex not in color_map:
        color_map[vertex] = 1

# Print the color map
print(color_map)