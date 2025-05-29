# Define the graph properties
vertices = list(range(10))
edges = [(0, 5), (1, 7)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Function to find a valid color for a vertex
def find_color(vertex, edges, color_map, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for v, adj in edges:
        if v == vertex and adj in color_map:
            adjacent_colors.add(color_map[adj])
        if adj == vertex and v in color_map:
            adjacent_colors.add(color_map[v])
    
    # Assign the smallest available color
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_color(vertex, edges, color_map, colors)

# Print the color map
print(color_map)