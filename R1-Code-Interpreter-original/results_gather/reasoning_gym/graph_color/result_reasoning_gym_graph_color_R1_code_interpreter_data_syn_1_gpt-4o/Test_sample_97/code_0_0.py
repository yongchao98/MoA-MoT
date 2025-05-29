# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (2, 9)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Assign colors to vertices connected by edges
color_map[0] = 1
color_map[1] = 2
color_map[2] = 1
color_map[9] = 2

# Assign colors to remaining vertices
for vertex in vertices:
    if vertex not in color_map:
        color_map[vertex] = 1  # Assign any available color

print(color_map)