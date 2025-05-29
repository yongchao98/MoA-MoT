import json

def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring():
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    colors = [1, 2, 3]
    edges = [(0, 9), (2, 4), (3, 9), (5, 6), (6, 7)]
    
    # Create adjacency list
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize coloring with first available color
    coloring = {v: 1 for v in vertices}
    
    # Color vertices with constraints
    for vertex in vertices:
        # Get colors of neighbors
        neighbor_colors = set(coloring[n] for n in graph[vertex] if n in coloring)
        
        # Find first available color that's not used by neighbors
        for color in colors:
            if color not in neighbor_colors:
                coloring[vertex] = color
                break
    
    # Convert to string keys for JSON
    result = {str(k): v for k, v in coloring.items()}
    print(json.dumps(result))

graph_coloring()