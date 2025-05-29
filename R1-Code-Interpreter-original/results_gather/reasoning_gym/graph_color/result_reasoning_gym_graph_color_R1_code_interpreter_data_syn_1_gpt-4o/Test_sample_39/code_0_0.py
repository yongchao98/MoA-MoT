def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {vertex: None for vertex in vertices}
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of adjacent vertices
        adjacent_colors = set(color_assignment[adj] for adj in vertices if (vertex, adj) in edges or (adj, vertex) in edges)
        
        # Assign the smallest possible color that is not used by adjacent vertices
        for color in colors:
            if color not in adjacent_colors:
                color_assignment[vertex] = color
                break
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 7), (1, 9), (3, 4), (3, 6), (5, 6), (5, 7), (5, 8)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result
print(color_assignment)