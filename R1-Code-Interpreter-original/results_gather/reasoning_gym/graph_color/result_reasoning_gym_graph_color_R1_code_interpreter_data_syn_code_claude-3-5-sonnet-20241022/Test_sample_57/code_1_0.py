import json

def is_safe(graph, vertex, colors, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring():
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 5), (1, 9), (2, 5), (3, 5), (4, 5)]
    num_colors = 3
    
    # Create adjacency list representation
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    colors = {}
    for vertex in vertices:
        colors[vertex] = None
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in range(1, num_colors + 1):
            if is_safe(adj_list, vertex, colors, color):
                colors[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                colors[vertex] = None
        return False
    
    # Find coloring
    if backtrack(0):
        # Convert colors to strings for JSON
        result = {str(k): v for k, v in colors.items()}
        print(json.dumps(result))
    else:
        print("No solution exists")

# Solve the problem
graph_coloring()