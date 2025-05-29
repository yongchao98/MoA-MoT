def graph_coloring(vertices, edges, colors):
    # Initialize the color map with None
    color_map = {v: None for v in vertices}
    
    # Function to check if a color can be assigned to a vertex
    def can_color(vertex, color):
        for edge in edges:
            if vertex in edge:
                adjacent = edge[0] if edge[1] == vertex else edge[1]
                if color_map[adjacent] == color:
                    return False
        return True
    
    # Assign colors to each vertex
    for vertex in vertices:
        for color in colors:
            if can_color(vertex, color):
                color_map[vertex] = color
                break
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (1, 6), (2, 6), (4, 7), (7, 8)]
colors = [1, 2, 3]

# Get the coloring result
coloring_result = graph_coloring(vertices, edges, colors)

# Print the result
print(coloring_result)