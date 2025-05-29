# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 7), (5, 7)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Function to find the smallest available color for a vertex
def find_color(vertex, color_map, edges, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = {color_map.get(neighbor) for u, v in edges for neighbor in (u, v) if neighbor == vertex}
    # Return the first available color that is not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Color the vertices
for vertex in vertices:
    color_map[vertex] = find_color(vertex, color_map, edges, colors)

# Print the resulting color map
print(color_map)