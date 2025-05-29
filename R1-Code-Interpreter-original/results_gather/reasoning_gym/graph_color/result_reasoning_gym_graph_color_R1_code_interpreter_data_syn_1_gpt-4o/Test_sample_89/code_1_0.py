# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 2), (3, 9), (6, 8)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Function to find a valid color for a vertex
def find_color(vertex, color_map, edges):
    # Get the colors of adjacent vertices
    adjacent_colors = set()
    for u, v in edges:
        if u == vertex and v in color_map:
            adjacent_colors.add(color_map[v])
        elif v == vertex and u in color_map:
            adjacent_colors.add(color_map[u])
    # Find the first available color
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_color(vertex, color_map, edges)

# Print the color map
print(color_map)