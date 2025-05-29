# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (1, 7), (2, 8), (5, 9)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Function to get adjacent vertices
def get_adjacent(vertex, edges):
    return [v for e in edges for v in e if vertex in e and v != vertex]

# Assign colors to each vertex
for vertex in vertices:
    # Get the colors of adjacent vertices
    adjacent_colors = {color_map[v] for v in get_adjacent(vertex, edges) if v in color_map}
    # Find the smallest available color
    for color in colors:
        if color not in adjacent_colors:
            color_map[vertex] = color
            break

# Print the color map
print(color_map)