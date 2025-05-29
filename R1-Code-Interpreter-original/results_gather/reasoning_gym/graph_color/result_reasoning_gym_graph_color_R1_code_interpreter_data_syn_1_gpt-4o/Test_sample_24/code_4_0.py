def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment for each vertex
    color_assignment = {v: None for v in vertices}
    
    # Assign colors to each vertex
    for vertex in vertices:
        # Get the colors of adjacent vertices
        adjacent_colors = set()
        for edge in edges:
            if edge[0] == vertex and color_assignment[edge[1]] is not None:
                adjacent_colors.add(color_assignment[edge[1]])
            elif edge[1] == vertex and color_assignment[edge[0]] is not None:
                adjacent_colors.add(color_assignment[edge[0]])
        
        # Assign the smallest available color
        for color in colors:
            if color not in adjacent_colors:
                color_assignment[vertex] = color
                break
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 8), (1, 9), (3, 6), (4, 8), (5, 7), (5, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result as a JSON map
import json
print(json.dumps(color_assignment))