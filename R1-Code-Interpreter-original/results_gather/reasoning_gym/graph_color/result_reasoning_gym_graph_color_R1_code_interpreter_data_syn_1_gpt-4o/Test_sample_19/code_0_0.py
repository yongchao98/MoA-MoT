def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {v: None for v in vertices}
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of adjacent vertices
        adjacent_colors = set()
        for edge in edges:
            if edge[0] == vertex:
                adjacent_colors.add(color_assignment[edge[1]])
            elif edge[1] == vertex:
                adjacent_colors.add(color_assignment[edge[0]])
        
        # Assign the smallest possible color that is not used by adjacent vertices
        for color in colors:
            if color not in adjacent_colors:
                color_assignment[vertex] = color
                break
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 9), (2, 5), (3, 9), (4, 7), (4, 8)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result
print(color_assignment)