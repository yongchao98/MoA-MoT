def graph_coloring(vertices, edges, colors):
    # Initialize the color map with None
    color_map = {v: None for v in vertices}
    
    # Function to find a valid color for a vertex
    def find_color(vertex):
        # Get the colors of adjacent vertices
        adjacent_colors = {color_map[v] for u, v in edges if u == vertex and color_map[v] is not None}
        adjacent_colors.update({color_map[u] for u, v in edges if v == vertex and color_map[u] is not None})
        
        # Find the first available color
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
edges = [(0, 2), (0, 3), (2, 5), (3, 5), (4, 9)]
colors = [1, 2, 3]

# Get the coloring result
coloring_result = graph_coloring(vertices, edges, colors)

# Print the result as a JSON map
import json
print(json.dumps(coloring_result))