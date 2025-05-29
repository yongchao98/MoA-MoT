def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
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
    
    # Initialize graph as adjacency list
    graph_adj = {v: [] for v in vertices}
    for v1, v2 in [(1, 4), (1, 7), (2, 5), (2, 7), (3, 6)]:
        graph_adj[v1].append(v2)
        graph_adj[v2].append(v1)
    
    # Initialize coloring dictionary
    coloring = {}
    
    # Start backtracking from vertex 0
    backtrack(0)
    
    # Convert all values to integers in the result
    result = {str(k): v for k, v in coloring.items()}
    
    import json
    print(json.dumps(result))

# Define inputs
vertices = list(range(10))
colors = [1, 2, 3]

# Run the algorithm
graph_coloring(vertices, colors, vertices)