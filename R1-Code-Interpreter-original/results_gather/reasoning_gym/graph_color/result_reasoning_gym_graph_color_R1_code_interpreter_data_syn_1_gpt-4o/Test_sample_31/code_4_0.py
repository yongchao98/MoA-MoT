def is_safe(vertex, color, color_map, edges):
    for v, adj in edges:
        if v == vertex and color_map[adj] == color:
            return False
        if adj == vertex and color_map[v] == color:
            return False
    return True

def graph_coloring_util(vertices, edges, colors, color_map, vertex_index):
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    
    for color in colors:
        if is_safe(vertex, color, color_map, edges):
            color_map[vertex] = color
            if graph_coloring_util(vertices, edges, colors, color_map, vertex_index + 1):
                return True
            color_map[vertex] = None
    
    return False

def graph_coloring(vertices, edges, colors):
    color_map = {v: None for v in vertices}
    if graph_coloring_util(vertices, edges, colors, color_map, 0):
        return color_map
    else:
        return None

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (1, 4), (2, 3), (2, 4), (3, 8), (4, 8), (5, 7), (6, 7), (6, 8), (6, 9)]
colors = [1, 2, 3]

# Get the coloring solution
coloring_solution = graph_coloring(vertices, edges, colors)

# Print the solution
print(coloring_solution)