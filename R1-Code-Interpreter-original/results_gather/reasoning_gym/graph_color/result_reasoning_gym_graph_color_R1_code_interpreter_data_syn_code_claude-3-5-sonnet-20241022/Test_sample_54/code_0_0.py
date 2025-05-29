def solve_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(1, 6), (2, 8), (4, 6)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring (all vertices can start with color 1)
    coloring = {v: 1 for v in vertices}
    
    # For each vertex that's part of an edge, ensure proper coloring
    for v1, v2 in edges:
        # If adjacent vertices have same color, change one of them
        if coloring[v1] == coloring[v2]:
            # Try to find a different color for v2
            for color in colors:
                if color != coloring[v1]:
                    coloring[v2] = color
                    break
    
    # Convert to JSON string format
    result = "<<<" + str(coloring).replace("'", '"') + ">>>"
    print(result)

solve_graph_coloring()