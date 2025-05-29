def is_valid_color(graph, colors, vertex, color):
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, color_list):
    n = len(graph)
    colors = {}
    
    def backtrack(vertex):
        if vertex == n:
            return True
        
        for color in color_list:
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex + 1):
                    return True
                colors.pop(vertex, None)
        
        return False
    
    # Create adjacency list representation
    adj_list = {i: [] for i in range(n)}
    edges = [(0, 1), (0, 7), (0, 9), (1, 4), (2, 4), (3, 5), (3, 6), (6, 8), (7, 9)]
    for u, v in edges:
        adj_list[u].append(v)
        adj_list[v].append(u)
    
    # Find solution
    backtrack(0)
    
    # Convert to required format
    result = {str(k): v for k, v in colors.items()}
    print(f"<<<{result}>>>")

# Run the solution
vertices = list(range(10))
colors = [1, 2, 3]
graph_coloring(vertices, colors)