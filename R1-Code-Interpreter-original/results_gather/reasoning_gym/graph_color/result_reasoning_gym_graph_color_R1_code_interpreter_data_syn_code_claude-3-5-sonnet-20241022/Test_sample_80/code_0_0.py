def is_valid_color(graph, colors, vertex, color):
    # Check if any neighbor has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(vertices, edges, available_colors):
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    colors = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in available_colors:
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                colors.pop(vertex, None)
        
        return False
    
    # Start backtracking from first vertex
    backtrack(0)
    
    # Convert to string keys for JSON compatibility
    result = {str(k): v for k, v in colors.items()}
    import json
    print(json.dumps(result))

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 9), (2, 7), (2, 8), (3, 6), (3, 9), (7, 9)]
colors = [1, 2, 3]

graph_coloring(vertices, edges, colors)