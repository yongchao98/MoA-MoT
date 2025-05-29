def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {vertex: None for vertex in vertices}
    
    # Create an adjacency list from the edges
    adjacency_list = {vertex: [] for vertex in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    # Assign colors to each vertex
    for vertex in vertices:
        # Find the colors used by adjacent vertices
        used_colors = {color_assignment[neighbor] for neighbor in adjacency_list[vertex] if color_assignment[neighbor] is not None}
        
        # Assign the smallest available color
        for color in colors:
            if color not in used_colors:
                color_assignment[vertex] = color
                break
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 2), (1, 8), (2, 3), (3, 4), (8, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)
print(color_assignment)