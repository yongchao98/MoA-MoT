def graph_coloring(vertices, edges, colors):
    # Initialize the color map
    color_map = {v: None for v in vertices}
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of adjacent vertices
        adjacent_colors = set()
        for edge in edges:
            if edge[0] == vertex:
                if color_map[edge[1]] is not None:
                    adjacent_colors.add(color_map[edge[1]])
            elif edge[1] == vertex:
                if color_map[edge[0]] is not None:
                    adjacent_colors.add(color_map[edge[0]])
        
        # Assign the smallest possible color
        for color in colors:
            if color not in adjacent_colors:
                color_map[vertex] = color
                break
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (1, 8), (3, 6), (3, 9), (4, 7), (5, 9)]
colors = [1, 2, 3]

# Get the coloring
coloring = graph_coloring(vertices, edges, colors)

# Print the result
print(coloring)