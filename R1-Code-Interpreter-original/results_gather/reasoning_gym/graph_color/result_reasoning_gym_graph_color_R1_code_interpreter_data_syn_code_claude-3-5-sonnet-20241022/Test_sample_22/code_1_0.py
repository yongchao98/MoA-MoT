def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    if vertex not in graph:
        return True
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def find_coloring():
    # Initialize vertices and edges
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(1, 7), (2, 5), (3, 4), (3, 7), (4, 9)]
    colors = [1, 2, 3]
    
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize coloring with first color for all vertices
    coloring = {v: 1 for v in vertices}
    
    # For each vertex, try to find a valid color
    for vertex in vertices:
        # Try each color
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                break
    
    # Convert to string keys for JSON format
    result = {str(k): v for k, v in coloring.items()}
    
    import json
    print(json.dumps(result))

# Run the coloring algorithm
find_coloring()