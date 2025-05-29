# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 7), (2, 5), (3, 4), (3, 7), (4, 9)]
colors = [1, 2, 3]

# Initialize the color assignment dictionary
color_assignment = {vertex: None for vertex in vertices}

# Function to find the lowest available color for a vertex
def find_color(vertex, edges, color_assignment, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if vertex in edge:
            adjacent_vertex = edge[0] if edge[1] == vertex else edge[1]
            if color_assignment[adjacent_vertex] is not None:
                adjacent_colors.add(color_assignment[adjacent_vertex])
    
    # Assign the lowest available color
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_assignment[vertex] = find_color(vertex, edges, color_assignment, colors)

# Print the color assignment
print(color_assignment)