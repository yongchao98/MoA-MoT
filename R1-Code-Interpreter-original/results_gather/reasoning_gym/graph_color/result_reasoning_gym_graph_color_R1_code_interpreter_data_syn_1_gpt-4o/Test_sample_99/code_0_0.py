# Define the graph properties
vertices = list(range(10))
edges = [(0, 4), (3, 4), (7, 8)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Function to find the smallest available color for a vertex
def find_color(vertex, edges, color_map, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = {color_map.get(v) for u, v in edges if u == vertex or v == vertex}
    # Return the first available color
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_color(vertex, edges, color_map, colors)

# Print the resulting color map
print(color_map)