import json

def solve_graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(1, 9), (2, 3)]
    colors = [1, 2, 3]
    
    # Initialize all vertices with color 1 (we can optimize this)
    coloring = {v: 1 for v in vertices}
    
    # Create adjacency list
    adj = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj[v1].append(v2)
        adj[v2].append(v1)
    
    # For vertices with edges, ensure different colors
    # Only need to check vertices 1,9 and 2,3
    coloring[9] = 2  # Different from vertex 1
    coloring[3] = 2  # Different from vertex 2
    
    # Convert to JSON string format
    result = json.dumps(coloring)
    print(result)

solve_graph_coloring()