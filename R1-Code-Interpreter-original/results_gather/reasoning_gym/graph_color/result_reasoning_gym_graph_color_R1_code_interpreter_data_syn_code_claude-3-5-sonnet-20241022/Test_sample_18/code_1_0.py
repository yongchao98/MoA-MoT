def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    if vertex not in graph:
        return True
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def find_coloring(vertices, edges, colors):
    # Create adjacency list representation of the graph
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize coloring with the first color
    coloring = {v: colors[0] for v in vertices}
    
    # Try to find a valid coloring
    changed = True
    while changed:
        changed = False
        for v in vertices:
            # If current color conflicts with neighbors
            if not is_safe(graph, v, coloring[v], coloring):
                # Try each color
                for c in colors:
                    coloring[v] = c
                    if is_safe(graph, v, c, coloring):
                        changed = True
                        break
    
    # Convert to required format
    result = {str(k): v for k, v in coloring.items()}
    import json
    print(json.dumps(result))

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (1, 7), (2, 5), (2, 7), (3, 6)]
colors = [1, 2, 3]

# Find coloring
find_coloring(vertices, edges, colors)