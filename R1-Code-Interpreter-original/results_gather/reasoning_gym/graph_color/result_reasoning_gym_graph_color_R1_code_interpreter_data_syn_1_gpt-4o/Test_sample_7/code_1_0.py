# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 2), (0, 5), (1, 2), (1, 9), (2, 6), (5, 9)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Function to check if a color can be assigned to a vertex
def is_valid_color(vertex, color, color_map, edges):
    for edge in edges:
        if edge[0] == vertex and color_map[edge[1]] == color:
            return False
        if edge[1] == vertex and color_map[edge[0]] == color:
            return False
    return True

# Function to solve the coloring problem using backtracking
def color_graph(vertex_index, vertices, color_map, edges, colors):
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    for color in colors:
        if is_valid_color(vertex, color, color_map, edges):
            color_map[vertex] = color
            if color_graph(vertex_index + 1, vertices, color_map, edges, colors):
                return True
            color_map[vertex] = None
    
    return False

# Solve the coloring problem
color_graph(0, vertices, color_map, edges, colors)

# Print the color map
print(color_map)