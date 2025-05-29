# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (1, 6), (4, 7), (7, 8)]
possible_colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Create an adjacency list for the graph
adjacency_list = {vertex: [] for vertex in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Function to find the smallest available color for a vertex
def find_smallest_color(vertex):
    # Get the colors of adjacent vertices
    adjacent_colors = {color_map[adj] for adj in adjacency_list[vertex] if color_map[adj] is not None}
    # Find the smallest color not used by adjacent vertices
    for color in possible_colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_smallest_color(vertex)

# Print the resulting color map
print(color_map)