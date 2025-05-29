import json

def color_graph():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(3, 7), (4, 8)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {}
    
    # For connected vertices, ensure different colors
    # Color vertex 3 with color 1, vertex 7 with color 2
    coloring[3] = 1
    coloring[7] = 2
    
    # Color vertex 4 with color 1, vertex 8 with color 2
    coloring[4] = 1
    coloring[8] = 2
    
    # Color remaining vertices with color 1 (could be any color)
    for v in vertices:
        if v not in coloring:
            coloring[v] = 1
    
    # Convert to string format
    result = {str(k): v for k, v in coloring.items()}
    print(json.dumps(result))

color_graph()