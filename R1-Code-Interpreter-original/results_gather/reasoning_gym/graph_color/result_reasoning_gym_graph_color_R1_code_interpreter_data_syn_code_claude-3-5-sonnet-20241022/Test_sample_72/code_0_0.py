def solve_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 7), (5, 7)]
    colors = [1, 2, 3]
    
    # Initialize solution with all vertices colored as 1
    # (since most vertices are isolated)
    solution = {str(v): 1 for v in vertices}
    
    # Handle vertex 7 and its neighbors
    # Vertex 7 is connected to 0 and 5
    # So we'll color:
    # - vertices 0 and 5 with color 1
    # - vertex 7 with color 2
    solution['7'] = 2
    solution['0'] = 1
    solution['5'] = 1
    
    # Convert to JSON format
    result = str(solution).replace("'", '"')
    print(result)

solve_graph_coloring()