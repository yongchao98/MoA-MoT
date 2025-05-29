def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def color_graph(vertices, edges, available_colors):
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize colors dictionary
    colors = {}
    
    # For each vertex, try each color
    for vertex in vertices:
        for color in available_colors:
            if is_valid_coloring(graph, colors, vertex, color):
                colors[vertex] = color
                break
    
    # Convert to string format for JSON
    result = {str(k): v for k, v in colors.items()}
    print(result)

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 5), (2, 5), (5, 9), (8, 9)]
available_colors = [1, 2, 3]

color_graph(vertices, edges, available_colors)