def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def color_graph():
    # Initialize graph as adjacency list
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(1, 3), (1, 5), (5, 8), (7, 8)]
    possible_colors = [1, 2, 3]
    
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize colors dictionary
    colors = {}
    
    # For each vertex
    for vertex in vertices:
        # Try each color
        for color in possible_colors:
            if is_valid_coloring(graph, colors, vertex, color):
                colors[vertex] = color
                break
    
    # Convert to required format
    result = {str(k): v for k, v in colors.items()}
    print(f"<<<{result}>>>")

# Run the solution
color_graph()