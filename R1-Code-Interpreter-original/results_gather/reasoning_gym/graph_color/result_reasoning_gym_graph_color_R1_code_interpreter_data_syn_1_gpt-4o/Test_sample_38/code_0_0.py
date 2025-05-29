# Define the graph
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (1, 9), (2, 3), (2, 6)]
colors = [1, 2, 3]

# Initialize the color assignment dictionary
color_assignment = {vertex: None for vertex in vertices}

# Create an adjacency list for the graph
adjacency_list = {vertex: [] for vertex in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Function to find the smallest available color for a vertex
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

print(color_assignment)