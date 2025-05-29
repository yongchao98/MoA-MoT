# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 9), (2, 3), (2, 6), (3, 8), (5, 7)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Function to find a valid color for a vertex
def find_color(vertex, color_map, edges, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if vertex in edge:
            adjacent_vertex = edge[0] if edge[1] == vertex else edge[1]
            if color_map[adjacent_vertex] is not None:
                adjacent_colors.add(color_map[adjacent_vertex])
    
    # Assign the smallest possible color that is not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_color(vertex, color_map, edges, colors)

# Print the resulting color map
print(color_map)