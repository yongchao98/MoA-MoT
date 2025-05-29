# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(5, 7), (5, 9)]

# Initialize the color map
color_map = {}

# Assign colors to the connected vertices
color_map[5] = 1
color_map[7] = 2
color_map[9] = 3

# Assign color 1 to the remaining vertices
for vertex in vertices:
    if vertex not in color_map:
        color_map[vertex] = 1

# Print the color map
print(color_map)