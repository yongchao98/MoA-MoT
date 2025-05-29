def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {vertex: None for vertex in vertices}
    
    # Function to find the smallest available color for a vertex
    def find_color(vertex):
        # Get the colors of adjacent vertices
        adjacent_colors = {color_assignment[adj] for adj in vertices if (vertex, adj) in edges or (adj, vertex) in edges}
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
edges = [(0, 5), (2, 5), (5, 9), (8, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)

# Print the result
print(color_assignment)