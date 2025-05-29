def find_graph_coloring():
    # Initialize vertices and edges
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 6), (1, 2)]
    colors = [1, 2, 3]
    
    # Initialize coloring dictionary
    coloring = {}
    
    # Create adjacency list for easier checking
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Color each vertex
    for vertex in vertices:
        # Get colors of adjacent vertices
        adjacent_colors = {coloring[adj] for adj in adj_list[vertex] if adj in coloring}
        
        # Choose first available color
        for color in colors:
            if color not in adjacent_colors:
                coloring[vertex] = color
                break
    
    # Convert to JSON-like string format
    result = "<<<{"
    result += ", ".join(f'"{v}": {coloring[v]}' for v in sorted(vertices))
    result += "}>>>"
    
    print(result)

find_graph_coloring()