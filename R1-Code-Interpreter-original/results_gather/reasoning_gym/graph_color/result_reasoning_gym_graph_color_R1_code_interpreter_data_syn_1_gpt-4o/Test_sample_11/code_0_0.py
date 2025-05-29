def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {v: None for v in vertices}
    
    # Create an adjacency list from the edges
    adjacency_list = {v: [] for v in vertices}
    for (u, v) in edges:
        adjacency_list[u].append(v)
        adjacency_list[v].append(u)
    
    # Assign colors to each vertex
    for vertex in vertices:
        # Find the colors of the adjacent vertices
        adjacent_colors = {color_assignment[adj] for adj in adjacency_list[vertex] if color_assignment[adj] is not None}
        
        # Assign the smallest possible color
        for color in colors:
            if color not in adjacent_colors:
                color_assignment[vertex] = color
                break
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (1, 8), (3, 6), (3, 9), (4, 7), (5, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)
print(color_assignment)