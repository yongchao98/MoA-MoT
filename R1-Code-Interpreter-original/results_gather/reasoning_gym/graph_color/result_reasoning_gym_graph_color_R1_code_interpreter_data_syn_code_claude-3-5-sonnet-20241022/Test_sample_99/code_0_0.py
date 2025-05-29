def color_graph():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 4), (3, 4), (7, 8)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for u, v in edges:
        adj_list[u].append(v)
        adj_list[v].append(u)
    
    # Initialize coloring
    coloring = {v: 1 for v in vertices}  # Start with all vertices colored 1
    
    # Color vertices with constraints
    for v in vertices:
        neighbor_colors = {coloring[n] for n in adj_list[v]}
        # Choose the first available color not used by neighbors
        for color in colors:
            if color not in neighbor_colors:
                coloring[v] = color
                break
    
    # Convert to string format
    result = "<<<" + str(coloring).replace("'", '"') + ">>>"
    print(result)

color_graph()