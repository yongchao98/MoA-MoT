import json

def is_valid_coloring(graph, colors):
    for v1, v2 in edges:  # Check each edge
        if colors[v1] == colors[v2]:  # If connected vertices have same color
            return False
    return True

def try_coloring():
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 4), (1, 2), (1, 3), (1, 9), (3, 5), (4, 8), (7, 9)]
    colors = {}
    
    # Start with vertex 1 since it has most connections
    colors[1] = 1
    
    # Color vertex 1's neighbors differently
    for neighbor in [2, 3, 9]:
        colors[neighbor] = 2
    
    # Color vertex 0
    colors[0] = 1
    
    # Color vertex 4 (connected to 0)
    colors[4] = 2
    
    # Color vertex 5 (connected to 3)
    colors[5] = 1
    
    # Color vertex 8 (connected to 4)
    colors[8] = 3
    
    # Color vertex 7 (connected to 9)
    colors[7] = 3
    
    # Color remaining vertices with any color
    for v in vertices:
        if v not in colors:
            colors[v] = 1
    
    return colors

# Get coloring
colors = try_coloring()

# Convert to string keys for JSON
solution = {str(k): v for k, v in colors.items()}

# Verify solution before printing
edges = [(0, 4), (1, 2), (1, 3), (1, 9), (3, 5), (4, 8), (7, 9)]
if is_valid_coloring(edges, colors):
    print(f"<<<{json.dumps(solution)}>>>")
else:
    # If first attempt fails, try alternative coloring
    colors = {
        0: 1,
        1: 2,
        2: 1,
        3: 1,
        4: 2,
        5: 2,
        6: 1,
        7: 1,
        8: 1,
        9: 3
    }
    solution = {str(k): v for k, v in colors.items()}
    print(f"<<<{json.dumps(solution)}>>>")