def verify_coloring(coloring, edges):
    # Check each edge for conflicts
    for edge in edges:
        v1, v2 = edge
        if coloring[v1] == coloring[v2]:
            return False
    return True

def get_neighbors(vertex, edges):
    neighbors = set()
    for v1, v2 in edges:
        if v1 == vertex:
            neighbors.add(v2)
        elif v2 == vertex:
            neighbors.add(v1)
    return neighbors

def get_valid_colors(vertex, current_coloring, edges, all_colors):
    neighbors = get_neighbors(vertex, edges)
    used_colors = {current_coloring[n] for n in neighbors if n in current_coloring}
    return [c for c in all_colors if c not in used_colors]

def color_graph(vertices, edges, colors):
    # Sort vertices by degree (number of connections)
    vertex_degrees = {v: len(get_neighbors(v, edges)) for v in vertices}
    sorted_vertices = sorted(vertices, key=lambda v: vertex_degrees[v], reverse=True)
    
    coloring = {}
    
    # Color vertices one by one
    for vertex in sorted_vertices:
        valid_colors = get_valid_colors(vertex, coloring, edges, colors)
        if not valid_colors:  # If no valid colors available, start over
            return None
        coloring[vertex] = valid_colors[0]
    
    # Verify final coloring
    if verify_coloring(coloring, edges):
        return coloring
    return None

# Problem setup
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (0, 6), (2, 6), (2, 8), (3, 4), (3, 6), (5, 6), (5, 7)]
colors = [1, 2, 3]

# Find solution
solution = color_graph(vertices, edges, colors)

# Format and print solution
if solution:
    result = {str(k): v for k, v in solution.items()}
    print(result)
else:
    print("No valid coloring found")