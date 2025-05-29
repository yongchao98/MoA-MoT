def get_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(2, 3), (2, 8)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {}
    
    # For vertex 2 and its neighbors
    coloring[2] = 1  # Assign first color to vertex 2
    coloring[3] = 2  # Assign second color to vertex 3
    coloring[8] = 2  # Assign second color to vertex 8
    
    # For all other vertices, assign first color (they're not connected to anything)
    for v in vertices:
        if v not in coloring:
            coloring[v] = 1
    
    # Convert to JSON string format
    result = "{"
    for v in vertices:
        result += f"\"{v}\": {coloring[v]}"
        if v < vertices[-1]:
            result += ", "
    result += "}"
    
    print(f"<<<{result}>>>")

get_graph_coloring()