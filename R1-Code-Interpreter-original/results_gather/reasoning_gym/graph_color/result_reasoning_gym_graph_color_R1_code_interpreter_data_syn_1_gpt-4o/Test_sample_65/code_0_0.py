def graph_coloring(vertices, edges, colors):
    # Initialize the color mapping
    color_map = {vertex: None for vertex in vertices}
    
    # Create an adjacency list
    adjacency_list = {vertex: [] for vertex in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    # Assign colors to each vertex
    for vertex in vertices:
        # Find the colors of adjacent vertices
        adjacent_colors = {color_map[adj] for adj in adjacency_list[vertex] if color_map[adj] is not None}
        
        # Assign the smallest available color
        for color in colors:
            if color not in adjacent_colors:
                color_map[vertex] = color
                break
    
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 5), (0, 6), (0, 8), (1, 3), (3, 9), (4, 5), (4, 6), (6, 8), (7, 8)]
colors = [1, 2, 3]

# Get the color mapping
color_mapping = graph_coloring(vertices, edges, colors)

# Print the result
print(color_mapping)