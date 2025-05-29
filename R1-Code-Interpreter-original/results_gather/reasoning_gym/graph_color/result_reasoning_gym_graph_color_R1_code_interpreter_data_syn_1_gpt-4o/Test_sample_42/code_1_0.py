def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment for each vertex
    color_assignment = {v: None for v in vertices}
    
    # Create an adjacency list for the graph
    adjacency_list = {v: [] for v in vertices}
    for u, v in edges:
        adjacency_list[u].append(v)
        adjacency_list[v].append(u)
    
    # Function to check if the current color assignment is valid
    def is_valid(vertex, color):
        for neighbor in adjacency_list[vertex]:
            if color_assignment[neighbor] == color:
                return False
        return True
    
    # Function to assign colors using backtracking
    def assign_colors(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in colors:
            if is_valid(vertex, color):
                color_assignment[vertex] = color
                if assign_colors(vertex_index + 1):
                    return True
                color_assignment[vertex] = None
        
        return False
    
    # Start the coloring process
    assign_colors(0)
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 6), (3, 6), (3, 9), (4, 5)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)
print(color_assignment)