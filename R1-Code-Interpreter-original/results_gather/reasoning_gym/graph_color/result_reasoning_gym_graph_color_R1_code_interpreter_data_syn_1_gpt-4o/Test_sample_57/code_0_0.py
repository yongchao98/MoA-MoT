# Define the graph
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 5), (1, 9), (2, 5), (3, 5), (4, 5)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Function to check if a color can be assigned to a vertex
def can_color(vertex, color):
    for edge in edges:
        if vertex in edge:
            adjacent_vertex = edge[0] if edge[1] == vertex else edge[1]
            if color_map[adjacent_vertex] == color:
                return False
    return True

# Assign colors to each vertex
for vertex in vertices:
    for color in colors:
        if can_color(vertex, color):
            color_map[vertex] = color
            break

print(color_map)