def solve_graph_coloring():
    vertices = list(range(10))
    edges = [(1, 4), (4, 5), (5, 8)]
    colors = [1, 2, 3]
    
    # Initialize all vertices with color 1
    coloring = {v: 1 for v in vertices}
    
    # Create adjacency list
    adj = {v: [] for v in vertices}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    
    # Color connected vertices
    for v in vertices:
        if adj[v]:  # if vertex has any neighbors
            neighbor_colors = {coloring[u] for u in adj[v]}
            # Find first available color not used by neighbors
            for color in colors:
                if color not in neighbor_colors:
                    coloring[v] = color
                    break
    
    # Convert to JSON string format
    result = str(coloring).replace("'", '"')
    print(f"<<<{result}>>>")

solve_graph_coloring()