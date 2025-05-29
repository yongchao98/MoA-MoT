def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {vertex: None for vertex in vertices}
    
    # Create an adjacency list
    adjacency_list = {vertex: [] for vertex in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    # Iterate over each vertex
    for vertex in vertices:
        # Find the colors of the adjacent vertices
        adjacent_colors = set()
        for neighbor in adjacency_list[vertex]:
            if color_assignment[neighbor] is not None:
                adjacent_colors.add(color_assignment[neighbor])
        
        # Assign the smallest available color
        for color in colors:
            if color not in adjacent_colors:
                color_assignment[vertex] = color
                break
    
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 5), (0, 7), (2, 3), (2, 7), (3, 5), (4, 7)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result
print(color_assignment)