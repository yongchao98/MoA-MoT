import json

def solve_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(4, 5), (4, 8), (5, 7)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {}
    
    # For vertices with edges, ensure proper coloring
    for vertex in vertices:
        if adj_list[vertex]:  # If vertex has connections
            # Get colors of neighbors
            neighbor_colors = {coloring.get(n) for n in adj_list[vertex]}
            # Choose first available color
            for color in colors:
                if color not in neighbor_colors:
                    coloring[vertex] = color
                    break
        else:  # If vertex has no connections
            coloring[vertex] = 1  # Use first color for isolated vertices
    
    # Fill in any uncolored vertices
    for vertex in vertices:
        if vertex not in coloring:
            coloring[vertex] = 1
    
    # Convert to string format
    result = {str(k): v for k, v in coloring.items()}
    return json.dumps(result)

print(solve_graph_coloring())