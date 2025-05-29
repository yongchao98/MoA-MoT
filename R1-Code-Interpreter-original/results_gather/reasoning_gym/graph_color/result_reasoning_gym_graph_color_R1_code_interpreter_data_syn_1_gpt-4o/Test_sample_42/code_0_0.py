def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment for each vertex
    color_assignment = {v: None for v in vertices}
    
    # Function to find the smallest available color for a vertex
    def find_color(vertex):
        # Get the colors of adjacent vertices
        adjacent_colors = {color_assignment[v] for u, v in edges if u == vertex or v == vertex}
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
edges = [(1, 6), (3, 6), (3, 9), (4, 5)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)
print(color_assignment)