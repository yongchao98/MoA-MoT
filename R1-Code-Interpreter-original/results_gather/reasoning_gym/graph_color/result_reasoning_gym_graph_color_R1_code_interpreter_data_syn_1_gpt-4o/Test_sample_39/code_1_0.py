def is_safe(vertex, color, color_assignment, edges):
    for edge in edges:
        if edge[0] == vertex and color_assignment[edge[1]] == color:
            return False
        if edge[1] == vertex and color_assignment[edge[0]] == color:
            return False
    return True

def graph_coloring_util(vertices, edges, colors, color_assignment, vertex_index):
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    for color in colors:
        if is_safe(vertex, color, color_assignment, edges):
            color_assignment[vertex] = color
            if graph_coloring_util(vertices, edges, colors, color_assignment, vertex_index + 1):
                return True
            color_assignment[vertex] = None
    
    return False

def graph_coloring(vertices, edges, colors):
    color_assignment = {vertex: None for vertex in vertices}
    if graph_coloring_util(vertices, edges, colors, color_assignment, 0):
        return color_assignment
    else:
        return None

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 7), (1, 9), (3, 4), (3, 6), (5, 6), (5, 7), (5, 8)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result
print(color_assignment)