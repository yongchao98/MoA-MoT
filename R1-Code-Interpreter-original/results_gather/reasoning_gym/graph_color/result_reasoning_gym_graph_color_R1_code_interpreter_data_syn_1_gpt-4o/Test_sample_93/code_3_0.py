# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 8), (3, 9), (4, 7)]
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
    # Find a color that is not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to vertices connected by edges first
for u, v in edges:
    if u not in color_map:
        color_map[u] = find_color(u, color_map, edges)
    if v not in color_map:
        color_map[v] = find_color(v, color_map, edges)

# Assign colors to the remaining vertices
for vertex in vertices:
    if vertex not in color_map:
        color_map[vertex] = find_color(vertex, color_map, edges)

# Print the color map
print(color_map)