import json

def get_neighbors(vertex, edges):
    neighbors = set()
    for edge in edges:
        if edge[0] == vertex:
            neighbors.add(edge[1])
        elif edge[1] == vertex:
            neighbors.add(edge[0])
    return neighbors

def find_valid_coloring():
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 6), (0, 8), (1, 4), (1, 6), (1, 8), (5, 6), (5, 8)]
    colors = [1, 2, 3]
    
    # Create adjacency map
    adjacency = {v: get_neighbors(v, edges) for v in vertices}
    
    # Sort vertices by degree (most connected first)
    vertices_by_degree = sorted(vertices, key=lambda v: len(adjacency[v]), reverse=True)
    
    # Initialize coloring
    coloring = {}
    
    # Color vertices one by one
    for vertex in vertices_by_degree:
        # Get colors used by neighbors
        neighbor_colors = {coloring.get(n) for n in adjacency[vertex]}
        
        # Find first available color
        for color in colors:
            if color not in neighbor_colors:
                coloring[vertex] = color
                break
    
    # Fill in remaining vertices with color 1 (for those with no connections)
    for vertex in vertices:
        if vertex not in coloring:
            coloring[vertex] = 1
    
    # Verify solution
    for edge in edges:
        if coloring[edge[0]] == coloring[edge[1]]:
            raise Exception(f"Invalid coloring: vertices {edge[0]} and {edge[1]} have same color")
    
    # Convert to string keys for JSON
    result = {str(k): v for k, v in coloring.items()}
    return result

# Find and print solution
result = find_valid_coloring()
print(f"<<<{json.dumps(result)}>>>")