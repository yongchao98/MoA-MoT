def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {v: None for v in vertices}
    
    # Create an adjacency list for the graph
    adjacency_list = {v: [] for v in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    # Color the vertices involved in edges first
    for vertex in vertices:
        if adjacency_list[vertex]:  # If the vertex has neighbors
            # Find the colors of adjacent vertices
            adjacent_colors = set(color_assignment[adj] for adj in adjacency_list[vertex] if color_assignment[adj] is not None)
            
            # Assign the smallest possible color that is not used by adjacent vertices
            for color in colors:
                if color not in adjacent_colors:
                    color_assignment[vertex] = color
                    break
    
    # Color the remaining vertices
    for vertex in vertices:
        if color_assignment[vertex] is None:
            color_assignment[vertex] = colors[0]  # Assign the first color available
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 9), (2, 3), (3, 7), (5, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result
print(color_assignment)