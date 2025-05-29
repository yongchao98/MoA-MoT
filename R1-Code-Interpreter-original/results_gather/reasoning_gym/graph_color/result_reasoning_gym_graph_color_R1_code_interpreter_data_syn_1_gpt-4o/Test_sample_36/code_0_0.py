# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 7), (2, 4), (6, 7)]
colors = [1, 2, 3]

# Initialize the color map
color_map = {vertex: None for vertex in vertices}

# Create an adjacency list for the graph
adjacency_list = {vertex: [] for vertex in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Greedy coloring algorithm
for vertex in vertices:
    # Find the colors of the adjacent vertices
    adjacent_colors = {color_map[adj] for adj in adjacency_list[vertex] if color_map[adj] is not None}
    
    # Assign the smallest available color
    for color in colors:
        if color not in adjacent_colors:
            color_map[vertex] = color
            break

# Print the color map
print(color_map)