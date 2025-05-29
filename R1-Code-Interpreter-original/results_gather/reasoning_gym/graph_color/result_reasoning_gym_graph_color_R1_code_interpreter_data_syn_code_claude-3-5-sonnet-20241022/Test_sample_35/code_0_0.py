def get_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(6, 8), (8, 9)]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {v: 1 for v in vertices}  # Default all to color 1
    
    # Color the connected component (6-8-9)
    coloring[6] = 1
    coloring[8] = 2
    coloring[9] = 1
    
    # Convert to JSON string format
    result = "<<<" + str(coloring).replace("'", '"') + ">>>"
    print(result)

get_graph_coloring()