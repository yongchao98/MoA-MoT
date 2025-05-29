import json

def color_graph():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 4), (0, 8), (3, 6)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring
    coloring = {}
    
    # Color vertices
    for vertex in vertices:
        # Get colors of adjacent vertices
        adjacent_colors = {coloring.get(adj) for adj in adj_list[vertex]}
        # Choose first available color
        for color in colors:
            if color not in adjacent_colors:
                coloring[vertex] = color
                break
    
    # Convert to string format
    result = {str(k): v for k, v in coloring.items()}
    print(f"<<<{json.dumps(result)}>>>")

# Run the coloring algorithm
color_graph()