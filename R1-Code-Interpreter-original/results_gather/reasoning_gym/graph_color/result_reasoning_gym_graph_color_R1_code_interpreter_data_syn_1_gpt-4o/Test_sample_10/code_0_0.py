# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (0, 7), (0, 9), (1, 4), (2, 4), (3, 5), (3, 6), (6, 8), (7, 9)]
possible_colors = [1, 2, 3]

# Initialize the color map
color_map = {}

# Create an adjacency list for the graph
adjacency_list = {v: [] for v in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Function to find the smallest available color
def find_smallest_color(used_colors):
    for color in possible_colors:
        if color not in used_colors:
            return color
    return None

# Assign colors to each vertex
for vertex in vertices:
    # Find colors of adjacent vertices
    used_colors = {color_map[adj] for adj in adjacency_list[vertex] if adj in color_map}
    # Assign the smallest possible color
    color_map[vertex] = find_smallest_color(used_colors)

# Print the result as a JSON map
import json
print(json.dumps(color_map))