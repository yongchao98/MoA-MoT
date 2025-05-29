# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 9), (2, 3), (2, 6), (3, 8), (5, 7)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Function to check if the current color assignment is valid
def is_valid(vertex, color, color_map, edges):
    for edge in edges:
        if vertex in edge:
            adjacent_vertex = edge[0] if edge[1] == vertex else edge[1]
            if color_map[adjacent_vertex] == color:
                return False
    return True

# Backtracking function to color the graph
def color_graph(vertex_index, color_map, edges, colors):
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    for color in colors:
        if is_valid(vertex, color, color_map, edges):
            color_map[vertex] = color
            if color_graph(vertex_index + 1, color_map, edges, colors):
                return True
            color_map[vertex] = None
    
    return False

# Start coloring from the first vertex
color_graph(0, color_map, edges, colors)

# Print the resulting color map
print(color_map)