def graph_coloring(vertices, edges, colors):
    # Initialize the color map with None for each vertex
    color_map = {vertex: None for vertex in vertices}
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of adjacent vertices
        adjacent_colors = set()
        for edge in edges:
            if edge[0] == vertex:
                adjacent_colors.add(color_map[edge[1]])
            elif edge[1] == vertex:
                adjacent_colors.add(color_map[edge[0]])
        
        # Assign the smallest available color
        for color in colors:
            if color not in adjacent_colors:
                color_map[vertex] = color
                break
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (1, 4), (2, 3), (2, 4), (3, 8), (4, 8), (5, 7), (6, 7), (6, 8), (6, 9)]
colors = [1, 2, 3]

# Get the coloring solution
coloring_solution = graph_coloring(vertices, edges, colors)

# Print the solution
print(coloring_solution)