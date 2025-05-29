# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 9), (2, 3), (2, 6), (3, 8), (5, 7)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Manually assign colors to ensure no two adjacent vertices have the same color
color_map[0] = 1
color_map[9] = 2
color_map[2] = 1
color_map[3] = 2
color_map[6] = 2
color_map[8] = 3
color_map[5] = 1
color_map[7] = 3
color_map[1] = 2
color_map[4] = 3

# Print the resulting color map
print(color_map)