def graph_coloring(vertices, edges, colors):
    # Initialize the color mapping
    color_map = {v: None for v in vertices}
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of the adjacent vertices
        adjacent_colors = {color_map[adj] for adj in vertices if (vertex, adj) in edges or (adj, vertex) in edges}
        
        # Assign the smallest possible color
        for color in colors:
            if color not in adjacent_colors:
                color_map[vertex] = color
                break
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (1, 2), (1, 6), (4, 6), (5, 6)]
colors = [1, 2, 3]

# Get the color mapping
color_mapping = graph_coloring(vertices, edges, colors)

# Print the result
print(color_mapping)