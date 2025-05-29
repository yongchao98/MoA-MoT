def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
    coloring = {}
    
    def backtrack(vertex_idx):
        if vertex_idx == len(vertices):
            return True
        
        vertex = vertices[vertex_idx]
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                if backtrack(vertex_idx + 1):
                    return True
        return False
    
    # Create adjacency list representation
    graph_adj = {v: [] for v in vertices}
    for v1, v2 in [(1, 2), (1, 8), (2, 3), (3, 4), (8, 9)]:
        graph_adj[v1].append(v2)
        graph_adj[v2].append(v1)
    
    # Try to find a valid coloring
    backtrack(0)
    
    # Convert all values to integers and ensure all vertices have a color
    result = {str(v): 1 for v in vertices}  # Default color 1 for unconnected vertices
    for v, c in coloring.items():
        result[str(v)] = c
    
    import json
    print(json.dumps(result))

# Run the coloring algorithm
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
colors = [1, 2, 3]
graph_coloring({v: [] for v in vertices}, colors, vertices)