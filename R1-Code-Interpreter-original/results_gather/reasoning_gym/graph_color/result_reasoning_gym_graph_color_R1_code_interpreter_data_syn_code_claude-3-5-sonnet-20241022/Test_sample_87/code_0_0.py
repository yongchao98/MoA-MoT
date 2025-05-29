def get_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 7), (2, 7), (4, 9)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {}
    
    # For each vertex
    for v in vertices:
        # Get colors of adjacent vertices
        adj_colors = {coloring.get(adj) for adj in adj_list[v]}
        # Choose first available color
        for color in colors:
            if color not in adj_colors:
                coloring[v] = color
                break
    
    # Convert to string format
    result = "<<<" + str({str(k): v for k, v in coloring.items()}) + ">>>"
    print(result)

get_graph_coloring()