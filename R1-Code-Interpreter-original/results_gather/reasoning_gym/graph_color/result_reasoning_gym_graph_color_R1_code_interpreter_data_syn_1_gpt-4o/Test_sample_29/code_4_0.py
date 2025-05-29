def graph_coloring(vertices, edges, colors):
    # Initialize the color map with None
    color_map = {v: None for v in vertices}
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of the adjacent vertices
        adjacent_colors = set()
        for edge in edges:
            if edge[0] == vertex and color_map[edge[1]] is not None:
                adjacent_colors.add(color_map[edge[1]])
            elif edge[1] == vertex and color_map[edge[0]] is not None:
                adjacent_colors.add(color_map[edge[0]])
        
        # Assign the smallest possible color
        for color in colors:
            if color not in adjacent_colors:
                color_map[vertex] = color
                break
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (0, 8), (1, 4), (1, 6), (1, 8), (5, 6), (5, 8)]
colors = [1, 2, 3]

# Get the color mapping
color_mapping = graph_coloring(vertices, edges, colors)

# Print the result
print(color_mapping)