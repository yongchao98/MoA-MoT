def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment for each vertex
    color_assignment = {vertex: None for vertex in vertices}
    
    # Function to find the smallest available color for a vertex
    def find_color(vertex):
        # Get the colors of adjacent vertices
        adjacent_colors = set()
        for edge in edges:
            if vertex in edge:
                adj_vertex = edge[0] if edge[1] == vertex else edge[1]
                if color_assignment[adj_vertex] is not None:
                    adjacent_colors.add(color_assignment[adj_vertex])
        
        # Find the smallest color not used by adjacent vertices
        for color in colors:
            if color not in adjacent_colors:
                return color
        return None
    
    # Assign colors to each vertex
    for vertex in vertices:
        color_assignment[vertex] = find_color(vertex)
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 2), (1, 5), (3, 5), (3, 9), (4, 5), (5, 7), (6, 8), (7, 9), (8, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result as a JSON map
import json
print(json.dumps(color_assignment))