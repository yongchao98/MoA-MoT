# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (0, 6), (1, 2), (3, 9)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Create an adjacency list for the graph
adjacency_list = {v: [] for v in vertices}
for u, v in edges:
    adjacency_list[u].append(v)
    adjacency_list[v].append(u)

# Assign colors to each vertex
for vertex in vertices:
    # Find the colors of adjacent vertices
    adjacent_colors = {color_map.get(neighbor) for neighbor in adjacency_list[vertex] if neighbor in color_map}
    
    # Assign the smallest available color
    for color in colors:
        if color not in adjacent_colors:
            color_map[vertex] = color
            break

# Print the color map
print(color_map)