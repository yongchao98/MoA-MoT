import json

def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
    coloring = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
        return False
    
    # Create adjacency list representation
    graph_adj = {v: [] for v in vertices}
    edges = [(0, 9), (2, 4), (3, 9), (5, 6), (6, 7)]
    for v1, v2 in edges:
        graph_adj[v1].append(v2)
        graph_adj[v2].append(v1)
    
    # Try to find a valid coloring
    backtrack(0)
    
    # Convert all values to integers for JSON
    result = {str(k): v for k, v in coloring.items()}
    print(json.dumps(result))

# Define the problem parameters
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
colors = [1, 2, 3]

# Solve the problem
graph_coloring({v: [] for v in vertices}, colors, vertices)