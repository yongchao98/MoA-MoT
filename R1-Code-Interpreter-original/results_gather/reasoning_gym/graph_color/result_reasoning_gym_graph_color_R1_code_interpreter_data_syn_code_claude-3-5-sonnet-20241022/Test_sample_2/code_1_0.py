import json

def solve_graph_coloring():
    # Create adjacency list representation
    adj = {i: [] for i in range(10)}
    edges = [(0, 9), (1, 8), (1, 9)]
    
    # Fill adjacency list
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    
    # Initialize colors
    colors = {i: 1 for i in range(10)}
    
    # Process vertices with edges first
    vertices_with_edges = set()
    for edge in edges:
        vertices_with_edges.add(edge[0])
        vertices_with_edges.add(edge[1])
    
    # Color vertices with edges
    for vertex in vertices_with_edges:
        neighbor_colors = {colors[n] for n in adj[vertex]}
        # Find first available color
        for color in [1, 2, 3]:
            if color not in neighbor_colors:
                colors[vertex] = color
                break
    
    # Convert to string keys for JSON format
    result = {str(k): v for k, v in colors.items()}
    
    # Verify solution
    for u, v in edges:
        if colors[u] == colors[v]:
            print(f"Invalid coloring found between vertices {u} and {v}!")
            return None
    
    print(json.dumps(result))

solve_graph_coloring()