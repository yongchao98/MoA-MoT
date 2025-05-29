# Define the graph
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 9), (2, 4), (3, 9), (5, 6), (6, 7)]
colors = [1, 2, 3]

# Initialize the color assignment
color_assignment = {vertex: None for vertex in vertices}

# Create an adjacency list
adjacency_list = {vertex: [] for vertex in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Function to find the smallest available color
def find_smallest_color(vertex):
    # Get the colors of adjacent vertices
    adjacent_colors = {color_assignment[adj] for adj in adjacency_list[vertex] if color_assignment[adj] is not None}
    # Find the smallest color not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_assignment[vertex] = find_smallest_color(vertex)

# Check for conflicts and adjust if necessary
for vertex in vertices:
    for neighbor in adjacency_list[vertex]:
        if color_assignment[vertex] == color_assignment[neighbor]:
            # Reassign a different color to the vertex
            color_assignment[vertex] = find_smallest_color(vertex)

print(color_assignment)