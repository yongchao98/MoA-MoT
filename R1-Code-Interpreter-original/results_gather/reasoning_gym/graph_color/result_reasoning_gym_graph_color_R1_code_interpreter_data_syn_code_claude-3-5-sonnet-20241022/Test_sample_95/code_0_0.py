def color_graph():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 4), (1, 4), (5, 9)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for u, v in edges:
        adj_list[u].append(v)
        adj_list[v].append(u)
    
    # Initialize coloring
    coloring = {}
    
    # Color vertices with constraints
    for v in vertices:
        # Get colors of adjacent vertices
        adj_colors = {coloring.get(adj) for adj in adj_list[v]}
        # Choose first available color
        for color in colors:
            if color not in adj_colors:
                coloring[v] = color
                break
    
    # Convert to required format
    result = {str(k): v for k, v in coloring.items()}
    print("<<<" + str(result) + ">>>")

color_graph()