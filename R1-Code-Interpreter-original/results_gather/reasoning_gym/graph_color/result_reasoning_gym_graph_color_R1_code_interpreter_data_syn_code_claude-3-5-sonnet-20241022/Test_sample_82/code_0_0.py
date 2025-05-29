import json

def color_graph():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 5), (1, 7)]
    colors = [1, 2, 3]
    
    # Initialize coloring dictionary
    coloring = {}
    
    # Create adjacency list for easier checking
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Color each vertex
    for vertex in vertices:
        # Get colors of adjacent vertices
        adjacent_colors = {coloring.get(adj) for adj in adj_list[vertex]}
        
        # Choose the first available color
        for color in colors:
            if color not in adjacent_colors:
                coloring[vertex] = color
                break
    
    # Convert to string format required
    result = {str(k): v for k, v in coloring.items()}
    return f"<<<{json.dumps(result)}>>>"

print(color_graph())