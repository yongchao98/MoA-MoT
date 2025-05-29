# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 2), (3, 9), (6, 8)]
colors = [1, 2, 3]

# Initialize a dictionary to store the color of each vertex
vertex_colors = {}

# Function to find the smallest available color for a vertex
def find_color(vertex, edges, vertex_colors, colors):
    # Find colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if edge[0] == vertex and edge[1] in vertex_colors:
            adjacent_colors.add(vertex_colors[edge[1]])
        elif edge[1] == vertex and edge[0] in vertex_colors:
            adjacent_colors.add(vertex_colors[edge[0]])
    
    # Assign the smallest color not used by adjacent vertices
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    vertex_colors[vertex] = find_color(vertex, edges, vertex_colors, colors)

# Print the resulting color mapping
print(vertex_colors)