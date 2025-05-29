def is_valid_coloring(graph, colors, vertex, color):
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def color_graph_backtrack(graph, vertex, colors, available_colors):
    # If all vertices are colored, we're done
    if vertex == len(graph):
        return True
    
    # Skip vertices that are already colored
    if vertex in colors:
        return color_graph_backtrack(graph, vertex + 1, colors, available_colors)
    
    # Try each color for current vertex
    for color in available_colors:
        if is_valid_coloring(graph, colors, vertex, color):
            colors[vertex] = color
            if color_graph_backtrack(graph, vertex + 1, colors, available_colors):
                return True
            colors.pop(vertex, None)
    
    return False

def solve_graph_coloring():
    # Initialize graph as adjacency list
    graph = {i: [] for i in range(10)}
    edges = [(0, 6), (0, 7), (1, 7), (4, 7), (5, 7), (7, 8)]
    
    # Build adjacency list
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    colors = {}
    available_colors = [1, 2, 3]
    
    # Start with vertex 7 (most connected)
    colors[7] = 1  # Try first color for vertex 7
    
    # Color remaining vertices using backtracking
    color_graph_backtrack(graph, 0, colors, available_colors)
    
    # Ensure all vertices have a color (assign 1 to any uncolored vertices)
    for vertex in range(10):
        if vertex not in colors:
            colors[vertex] = 1
    
    # Format output
    result = "{"
    for i in range(10):
        result += f'"{i}": {colors[i]}'
        if i < 9:
            result += ", "
    result += "}"
    print(f"<<<{result}>>>")

solve_graph_coloring()