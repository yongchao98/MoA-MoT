def is_valid_coloring(graph, colors, vertex, color):
    # Check if assigning 'color' to 'vertex' is valid
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def color_graph():
    # Initialize graph as adjacency list
    V = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    E = [(1, 6), (3, 6), (3, 9), (4, 5)]
    available_colors = [1, 2, 3]
    
    # Create adjacency list
    graph = {v: [] for v in V}
    for v1, v2 in E:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize colors dictionary
    colors = {}
    
    # For each vertex
    for vertex in V:
        # Try each color
        for color in available_colors:
            if is_valid_coloring(graph, colors, vertex, color):
                colors[vertex] = color
                break
    
    # Convert to required format
    result = {str(k): v for k, v in colors.items()}
    print(f"<<<{result}>>>")

# Run the solution
color_graph()