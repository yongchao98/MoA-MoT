# Define the vertices and edges of the graph
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (1, 8), (2, 8)]
colors = [1, 2, 3]

# Initialize a dictionary to store the color of each vertex
vertex_colors = {vertex: None for vertex in vertices}

# Function to find a valid color for a vertex
def find_color(vertex, edges, vertex_colors, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if edge[0] == vertex:
            adjacent_colors.add(vertex_colors[edge[1]])
        elif edge[1] == vertex:
            adjacent_colors.add(vertex_colors[edge[0]])
    
    # Assign the smallest available color
    for color in colors:
        if color not in adjacent_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    vertex_colors[vertex] = find_color(vertex, edges, vertex_colors, colors)

# Print the resulting coloring
print(vertex_colors)