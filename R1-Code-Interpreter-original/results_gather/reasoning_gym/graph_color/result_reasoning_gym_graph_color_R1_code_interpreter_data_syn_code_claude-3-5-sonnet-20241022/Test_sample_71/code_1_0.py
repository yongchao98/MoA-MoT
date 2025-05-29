def get_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(1, 7), (3, 7), (6, 9)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {}
    
    # For vertices with edges, ensure proper coloring
    # Start with vertex 7 as it has most connections
    coloring[7] = 1      # Assign first color to vertex 7
    coloring[1] = 2      # Must be different from vertex 7
    coloring[3] = 3      # Must be different from vertex 7
    coloring[6] = 1      # Can use any color
    coloring[9] = 2      # Different color from vertex 6
    
    # For all remaining vertices, assign first color (they're isolated)
    for v in vertices:
        if v not in coloring:
            coloring[v] = 1
    
    # Convert to string format and print
    result = "{"
    for v in vertices:
        result += f"\"{v}\": {coloring[v]}"
        if v != vertices[-1]:
            result += ", "
    result += "}"
    print(f"<<<{result}>>>")

get_graph_coloring()