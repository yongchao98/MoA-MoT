def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment for each vertex
    color_assignment = {v: None for v in vertices}
    
    # Create an adjacency list for the graph
    adjacency_list = {v: [] for v in vertices}
    for u, v in edges:
        adjacency_list[u].append(v)
        adjacency_list[v].append(u)
    
    # Assign colors to each vertex that is part of an edge
    for vertex in vertices:
        if adjacency_list[vertex]:  # Only consider vertices that have edges
            # Get the colors of adjacent vertices
            adjacent_colors = {color_assignment[adj] for adj in adjacency_list[vertex] if color_assignment[adj] is not None}
            
            # Find the smallest color not used by adjacent vertices
            for color in colors:
                if color not in adjacent_colors:
                    color_assignment[vertex] = color
                    break
        else:
            # Assign any color to vertices not part of any edge
            color_assignment[vertex] = colors[0]
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(2, 7), (2, 8), (3, 8), (3, 9), (4, 8), (5, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result as a JSON map
import json
print(json.dumps(color_assignment))