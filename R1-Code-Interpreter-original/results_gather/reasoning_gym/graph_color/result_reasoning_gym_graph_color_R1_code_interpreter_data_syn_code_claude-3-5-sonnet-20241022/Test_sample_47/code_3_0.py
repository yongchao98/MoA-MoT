def is_valid_coloring(graph, colors, vertex, color):
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def backtrack(graph, vertices_to_color, colors, available_colors):
    if not vertices_to_color:  # If no more vertices to color
        return True
    
    vertex = vertices_to_color[0]  # Get next vertex to color
    remaining_vertices = vertices_to_color[1:]
    
    for color in available_colors:
        if is_valid_coloring(graph, colors, vertex, color):
            colors[vertex] = color
            if backtrack(graph, remaining_vertices, colors, available_colors):
                return True
            colors.pop(vertex)
    
    return False

def solve_graph_coloring():
    # Initialize graph
    graph = {i: [] for i in range(10)}
    edges = [(0, 6), (0, 7), (1, 7), (4, 7), (5, 7), (7, 8)]
    
    # Build adjacency list
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Calculate vertex degrees
    vertex_degrees = {v: len(graph[v]) for v in range(10)}
    
    # Sort vertices by degree (most connected first)
    vertices_to_color = sorted(range(10), key=lambda x: vertex_degrees[x], reverse=True)
    
    colors = {}
    available_colors = [1, 2, 3]
    
    # Apply backtracking
    if backtrack(graph, vertices_to_color, colors, available_colors):
        # Fill in any uncolored vertices with color 1
        for v in range(10):
            if v not in colors:
                colors[v] = 1
        
        # Format output
        result = "{"
        for i in range(10):
            result += f'"{i}": {colors[i]}'
            if i < 9:
                result += ", "
        result += "}"
        print(f"<<<{result}>>>")

solve_graph_coloring()