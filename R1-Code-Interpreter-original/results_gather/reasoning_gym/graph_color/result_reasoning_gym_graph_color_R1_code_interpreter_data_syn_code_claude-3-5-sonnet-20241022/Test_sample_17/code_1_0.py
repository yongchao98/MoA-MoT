def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring():
    # Initialize vertices and colors
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    colors = [1, 2, 3]
    edges = [(1, 2), (1, 8), (2, 3), (3, 4), (8, 9)]
    
    # Create adjacency list
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize coloring with default color 1 for isolated vertices
    coloring = {v: 1 for v in vertices}
    
    # Color vertices that are part of edges
    connected_vertices = set()
    for v1, v2 in edges:
        connected_vertices.add(v1)
        connected_vertices.add(v2)
    
    # Process connected vertices
    for vertex in connected_vertices:
        # Find available color
        used_colors = {coloring[neighbor] for neighbor in graph[vertex]}
        for color in colors:
            if color not in used_colors:
                coloring[vertex] = color
                break
    
    # Convert to required format
    result = {str(v): coloring[v] for v in vertices}
    import json
    print(json.dumps(result))

# Run the coloring algorithm
graph_coloring()