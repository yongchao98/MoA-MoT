def solve_graph_coloring():
    # Initialize graph
    graph = {i: [] for i in range(10)}
    edges = [(0, 6), (0, 7), (1, 7), (4, 7), (5, 7), (7, 8)]
    
    # Build adjacency list
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    # Initialize all vertices with color 1
    colors = {i: 1 for i in range(10)}
    
    # Start with vertex 7 (most connected)
    colors[7] = 1  # Give vertex 7 color 1
    
    # Color its neighbors with color 2
    for neighbor in graph[7]:
        colors[neighbor] = 2
    
    # Color vertex 6 with color 3 (since it's connected to vertex 0 which has color 2)
    colors[6] = 3
    
    # Verify the coloring
    for v1, v2 in edges:
        if colors[v1] == colors[v2]:
            # If we find any conflict, adjust the coloring
            if v1 == 7 or v2 == 7:
                other = v2 if v1 == 7 else v1
                colors[other] = 2
            elif v1 == 6 or v2 == 6:
                other = v2 if v1 == 6 else v1
                colors[other] = 2
    
    # Format output
    result = "{"
    for i in range(10):
        result += f'"{i}": {colors[i]}'
        if i < 9:
            result += ", "
    result += "}"
    print(f"<<<{result}>>>")

solve_graph_coloring()