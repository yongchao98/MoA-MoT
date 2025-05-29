def graph_coloring(vertices, edges, colors):
    # Initialize the color map with None
    color_map = {v: None for v in vertices}
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of the adjacent vertices
        adjacent_colors = set()
        for v, adj in edges:
            if v == vertex and color_map[adj] is not None:
                adjacent_colors.add(color_map[adj])
            elif adj == vertex and color_map[v] is not None:
                adjacent_colors.add(color_map[v])
        
        # Assign the smallest possible color
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