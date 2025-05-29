# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 6), (2, 8), (4, 6)]
colors = [1, 2, 3]

# Initialize the color mapping
color_map = {}

# Function to find the smallest available color for a vertex
def find_color(vertex, edges, color_map, colors):
    # Find the colors of adjacent vertices
    adjacent_colors = set()
    for edge in edges:
        if vertex in edge:
            adjacent_vertex = edge[0] if edge[1] == vertex else edge[1]
            if adjacent_vertex in color_map:
                adjacent_colors.add(color_map[adjacent_vertex])
    
    # Assign the smallest available color
    for color in colors:
        if color not in adjacent_colors:
            return color

# Assign colors to each vertex
for vertex in vertices:
    color_map[vertex] = find_color(vertex, edges, color_map, colors)

# Print the resulting color map
print(color_map)