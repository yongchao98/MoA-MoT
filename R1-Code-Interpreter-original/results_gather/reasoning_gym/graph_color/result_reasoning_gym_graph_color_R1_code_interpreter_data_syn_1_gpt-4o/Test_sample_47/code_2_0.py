def is_safe(vertex, color, color_assignment, adjacency_list):
    for neighbor in adjacency_list[vertex]:
        if color_assignment[neighbor] == color:
            return False
    return True

def graph_coloring_util(vertices, adjacency_list, colors, color_assignment, vertex_index):
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    
    for color in colors:
        if is_safe(vertex, color, color_assignment, adjacency_list):
            color_assignment[vertex] = color
            if graph_coloring_util(vertices, adjacency_list, colors, color_assignment, vertex_index + 1):
                return True
            color_assignment[vertex] = None
    
    return False

def graph_coloring(vertices, edges, colors):
    color_assignment = {vertex: None for vertex in vertices}
    adjacency_list = {vertex: [] for vertex in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    if graph_coloring_util(vertices, adjacency_list, colors, color_assignment, 0):
        return color_assignment
    else:
        return None

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (0, 7), (1, 7), (4, 7), (5, 7), (7, 8)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result
print(color_assignment)