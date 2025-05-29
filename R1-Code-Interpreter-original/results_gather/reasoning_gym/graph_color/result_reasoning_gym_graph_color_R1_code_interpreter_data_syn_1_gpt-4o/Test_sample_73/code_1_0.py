def graph_coloring(vertices, edges, colors):
    # Initialize the color map
    color_map = {v: None for v in vertices}
    
    # Function to find the smallest available color for a vertex
    def find_color(vertex):
        # Get the colors of adjacent vertices
        adjacent_colors = set()
        for u, v in edges:
            if u == vertex and color_map[v] is not None:
                adjacent_colors.add(color_map[v])
            elif v == vertex and color_map[u] is not None:
                adjacent_colors.add(color_map[u])
        
        # Find the smallest color not used by adjacent vertices
        for color in colors:
            if color not in adjacent_colors:
                return color
        return None
    
    # Assign colors to each vertex
    for vertex in vertices:
        color_map[vertex] = find_color(vertex)
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 3), (1, 2), (3, 4)]
colors = [1, 2, 3]

# Get the coloring
coloring = graph_coloring(vertices, edges, colors)
print(coloring)