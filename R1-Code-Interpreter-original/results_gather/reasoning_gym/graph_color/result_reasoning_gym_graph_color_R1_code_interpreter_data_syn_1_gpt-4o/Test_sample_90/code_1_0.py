# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(2, 9), (5, 9)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Function to find the smallest available color for a vertex
def find_color(vertex, edges, color_map, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for u, v in edges:
        if u == vertex and v in color_map:
            adjacent_colors.add(color_map[v])
        elif v == vertex and u in color_map:
            adjacent_colors.add(color_map[u])
    # Return the first available color
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_color(vertex, edges, color_map, colors)

# Print the color map
print(color_map)