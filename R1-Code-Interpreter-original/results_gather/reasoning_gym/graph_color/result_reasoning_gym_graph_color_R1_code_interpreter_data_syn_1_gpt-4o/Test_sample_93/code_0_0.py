# Define the vertices and edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 8), (3, 9), (4, 7)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Create an adjacency list
adjacency_list = {vertex: [] for vertex in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Function to find the smallest available color
def find_smallest_color(used_colors):
    for color in colors:
        if color not in used_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    # Find the colors of adjacent vertices
    used_colors = {color_map[adj] for adj in adjacency_list[vertex] if color_map[adj] is not None}
    # Assign the smallest available color
    color_map[vertex] = find_smallest_color(used_colors)

# Print the color map
print(color_map)