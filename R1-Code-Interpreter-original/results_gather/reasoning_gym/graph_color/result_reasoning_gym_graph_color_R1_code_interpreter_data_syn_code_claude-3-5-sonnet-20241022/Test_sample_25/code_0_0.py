def get_graph_coloring():
    # Initialize all vertices with color 1 (any color would work for isolated vertices)
    coloring = {i: 1 for i in range(10)}
    
    # Handle the connected component (5, 7, 9)
    # Give vertex 5 color 1
    # Give vertices 7 and 9 color 2 (they can share a color since they're not connected)
    coloring[5] = 1
    coloring[7] = 2
    coloring[9] = 2
    
    # Convert to string format and return
    result = "{"
    for vertex in range(10):
        result += f'"{vertex}": {coloring[vertex]}'
        if vertex < 9:
            result += ", "
    result += "}"
    
    print(result)

get_graph_coloring()