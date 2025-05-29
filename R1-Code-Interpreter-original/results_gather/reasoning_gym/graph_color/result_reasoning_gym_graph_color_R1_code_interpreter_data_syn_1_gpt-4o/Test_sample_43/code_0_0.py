def graph_coloring(vertices, edges, colors):
    # Initialize the color map with None
    color_map = {v: None for v in vertices}
    
    # Function to check if the current color assignment is safe for vertex v
    def is_safe(v, c):
        for u in vertices:
            if (v, u) in edges or (u, v) in edges:
                if color_map[u] == c:
                    return False
        return True
    
    # Assign colors to vertices
    for v in vertices:
        for c in colors:
            if is_safe(v, c):
                color_map[v] = c
                break
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(4, 5), (4, 8), (5, 7)]
colors = [1, 2, 3]

# Get the color mapping
coloring = graph_coloring(vertices, edges, colors)
print(coloring)