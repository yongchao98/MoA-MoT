def is_safe(vertex, color, color_map, adjacency_list):
    # Check if the current color assignment is safe for the vertex
    for neighbor in adjacency_list[vertex]:
        if color_map[neighbor] == color:
            return False
    return True

def graph_coloring_backtrack(vertices, edges, colors):
    # Initialize the color mapping
    color_map = {vertex: None for vertex in vertices}
    
    # Create an adjacency list
    adjacency_list = {vertex: [] for vertex in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        
        for color in colors:
            if is_safe(vertex, color, color_map, adjacency_list):
                color_map[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                color_map[vertex] = None
        
        return False
    
    if backtrack(0):
        return color_map
    else:
        return None

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 5), (0, 6), (0, 8), (1, 3), (3, 9), (4, 5), (4, 6), (6, 8), (7, 8)]
colors = [1, 2, 3]

# Get the color mapping
color_mapping = graph_coloring_backtrack(vertices, edges, colors)

# Print the result
print(color_mapping)