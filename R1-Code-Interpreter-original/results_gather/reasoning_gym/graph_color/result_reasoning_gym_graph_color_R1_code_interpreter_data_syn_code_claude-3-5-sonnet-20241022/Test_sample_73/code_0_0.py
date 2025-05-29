def solve_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 3), (1, 2), (3, 4)]
    colors = [1, 2, 3]
    
    # Initialize all vertices with color 1
    coloring = {v: 1 for v in vertices}
    
    # Create adjacency list
    adj = {v: [] for v in vertices}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    
    # Color vertices that are part of edges
    for u, v in edges:
        # If adjacent vertices have same color
        if coloring[u] == coloring[v]:
            # Change color of second vertex
            coloring[v] = 2
    
    # Convert to required format
    result = {str(k): v for k, v in coloring.items()}
    print("<<<" + str(result) + ">>>")

solve_graph_coloring()