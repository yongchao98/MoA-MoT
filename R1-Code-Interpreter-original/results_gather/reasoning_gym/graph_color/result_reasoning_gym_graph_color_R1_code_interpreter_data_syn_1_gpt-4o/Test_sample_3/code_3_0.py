def is_valid(vertex, color, color_assignment, adjacency_list):
    for neighbor in adjacency_list[vertex]:
        if color_assignment[neighbor] == color:
            return False
    return True

def graph_coloring_backtrack(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {vertex: None for vertex in vertices}
    
    # Create an adjacency list from the edges
    adjacency_list = {vertex: [] for vertex in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    def assign_colors(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in colors:
            if is_valid(vertex, color, color_assignment, adjacency_list):
                color_assignment[vertex] = color
                if assign_colors(vertex_index + 1):
                    return True
                color_assignment[vertex] = None
        
        return False
    
    assign_colors(0)
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 9), (2, 5), (2, 7), (3, 6), (3, 8), (4, 5)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring_backtrack(vertices, edges, colors)
print(color_assignment)