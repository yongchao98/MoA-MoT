def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {vertex: None for vertex in vertices}
    
    # Function to check if a color can be assigned to a vertex
    def can_color(vertex, color):
        for edge in edges:
            if vertex in edge:
                adjacent_vertex = edge[0] if edge[1] == vertex else edge[1]
                if color_assignment[adjacent_vertex] == color:
                    return False
        return True
    
    # Assign colors to each vertex
    for vertex in vertices:
        for color in colors:
            if can_color(vertex, color):
                color_assignment[vertex] = color
                break
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (1, 2), (1, 3), (1, 9), (3, 5), (4, 8), (7, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)
print(color_assignment)