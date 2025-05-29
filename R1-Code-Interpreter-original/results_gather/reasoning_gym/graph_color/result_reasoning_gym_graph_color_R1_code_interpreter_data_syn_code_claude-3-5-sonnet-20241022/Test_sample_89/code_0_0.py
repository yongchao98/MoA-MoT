def solve_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(1, 2), (3, 9), (6, 8)]
    
    # Initialize all vertices with color 1 (since most are isolated)
    coloring = {v: 1 for v in vertices}
    
    # For each edge, ensure connected vertices have different colors
    for v1, v2 in edges:
        # If vertices have same color, change second vertex to color 2
        if coloring[v1] == coloring[v2]:
            coloring[v2] = 2
    
    # Convert to string format and return
    result = "<<<" + str(coloring).replace("'", '"') + ">>>"
    print(result)

solve_graph_coloring()