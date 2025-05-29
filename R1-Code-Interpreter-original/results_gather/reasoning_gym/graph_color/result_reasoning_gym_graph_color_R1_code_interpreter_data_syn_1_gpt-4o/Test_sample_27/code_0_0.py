# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (0, 6), (2, 6), (2, 8), (3, 4), (3, 6), (5, 6), (5, 7)]
colors = [1, 2, 3]

# Initialize the color assignment for each vertex
color_assignment = {vertex: 0 for vertex in vertices}

# Function to find the smallest available color for a vertex
def find_smallest_color(vertex, edges, color_assignment, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if edge[0] == vertex:
            adjacent_colors.add(color_assignment[edge[1]])
        elif edge[1] == vertex:
            adjacent_colors.add(color_assignment[edge[0]])
    
    # Find the smallest color not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_assignment[vertex] = find_smallest_color(vertex, edges, color_assignment, colors)

# Print the color assignment
print(color_assignment)