def is_valid_coloring(graph, colors, vertex, color):
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring():
    # Initialize graph as adjacency list
    graph = {i: [] for i in range(10)}
    edges = [(0, 6), (0, 7), (1, 7), (4, 7), (5, 7), (7, 8)]
    
    # Build adjacency list
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    vertices = list(range(10))
    available_colors = [1, 2, 3]
    colors = {}
    
    # For vertices not connected to any other vertex, assign color 1
    for vertex in vertices:
        if not graph[vertex]:
            colors[vertex] = 1
    
    # For connected vertices, try each color
    for vertex in vertices:
        if vertex not in colors:
            for color in available_colors:
                if is_valid_coloring(graph, colors, vertex, color):
                    colors[vertex] = color
                    break
    
    # Convert to JSON-like string
    result = "{"
    for i in range(10):
        result += f'"{i}": {colors[i]}'
        if i < 9:
            result += ", "
    result += "}"
    print(f"<<<{result}>>>")

graph_coloring()