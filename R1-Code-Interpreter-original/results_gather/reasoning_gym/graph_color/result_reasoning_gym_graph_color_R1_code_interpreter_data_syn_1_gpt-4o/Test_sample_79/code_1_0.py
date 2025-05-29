# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(2, 3), (2, 8)]
colors = [1, 2, 3]

# Initialize a dictionary to store the color of each vertex
vertex_colors = {}

# Function to find a valid color for a vertex
def find_color(vertex, edges, vertex_colors, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for v, adj in edges:
        if v == vertex and adj in vertex_colors:
            adjacent_colors.add(vertex_colors[adj])
        elif adj == vertex and v in vertex_colors:
            adjacent_colors.add(vertex_colors[v])
    # Return the first available color that is not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    vertex_colors[vertex] = find_color(vertex, edges, vertex_colors, colors)

# Print the resulting color mapping
print(vertex_colors)