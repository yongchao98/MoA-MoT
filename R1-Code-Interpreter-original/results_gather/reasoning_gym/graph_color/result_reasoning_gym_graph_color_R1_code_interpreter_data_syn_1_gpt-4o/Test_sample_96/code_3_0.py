# Define the graph
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (4, 5), (5, 8)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Function to find the smallest available color for a vertex
def find_color(vertex, edges, color_map, colors):
    # Find colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if edge[0] == vertex:
            if color_map[edge[1]] is not None:
                adjacent_colors.add(color_map[edge[1]])
        elif edge[1] == vertex:
            if color_map[edge[0]] is not None:
                adjacent_colors.add(color_map[edge[0]])
    
    # Find the smallest color not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_color(vertex, edges, color_map, colors)

# Print the color map
print(color_map)