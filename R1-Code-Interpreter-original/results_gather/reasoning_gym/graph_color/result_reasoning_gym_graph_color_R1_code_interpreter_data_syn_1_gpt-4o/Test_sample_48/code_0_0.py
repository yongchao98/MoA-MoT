# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 7), (1, 8), (4, 7), (5, 9), (8, 9)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Function to find the smallest available color for a vertex
def find_color(vertex, edges, color_map, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if vertex in edge:
            adjacent_vertex = edge[0] if edge[1] == vertex else edge[1]
            if color_map[adjacent_vertex] is not None:
                adjacent_colors.add(color_map[adjacent_vertex])
    
    # Find the smallest color not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Color the vertices
for vertex in vertices:
    color_map[vertex] = find_color(vertex, edges, color_map, colors)

print(color_map)