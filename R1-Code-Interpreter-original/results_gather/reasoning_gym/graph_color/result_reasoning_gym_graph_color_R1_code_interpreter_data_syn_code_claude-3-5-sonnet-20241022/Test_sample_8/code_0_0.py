import json

def solve_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(1, 5), (1, 9), (3, 6), (7, 9)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {}
    
    # For vertices with edges, ensure no conflicts
    for vertex in vertices:
        if adj_list[vertex]:  # If vertex has connections
            # Get colors of neighbors
            neighbor_colors = {coloring.get(n) for n in adj_list[vertex]}
            # Choose first available color
            for color in colors:
                if color not in neighbor_colors:
                    coloring[vertex] = color
                    break
        else:
            # If no connections, use first color
            coloring[vertex] = 1
    
    # Fill in any uncolored vertices
    for vertex in vertices:
        if vertex not in coloring:
            coloring[vertex] = 1
    
    # Convert to string format required
    result = {str(k): v for k, v in coloring.items()}
    print(f"<<<{json.dumps(result)}>>>")

solve_graph_coloring()