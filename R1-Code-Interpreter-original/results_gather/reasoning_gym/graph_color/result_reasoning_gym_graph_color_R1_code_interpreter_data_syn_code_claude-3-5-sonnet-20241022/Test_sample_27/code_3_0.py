def is_safe(vertex, color, coloring, edges):
    for v1, v2 in edges:
        if v1 == vertex and v2 in coloring and coloring[v2] == color:
            return False
        if v2 == vertex and v1 in coloring and coloring[v1] == color:
            return False
    return True

def solve_coloring(vertices, edges, colors):
    # Initialize coloring with vertex 6 (most connected)
    coloring = {6: 1}  # Start with vertex 6 as color 1
    uncolored = set(vertices) - {6}
    
    # Color vertices connected to 6 first
    connected_to_6 = set()
    for v1, v2 in edges:
        if v1 == 6:
            connected_to_6.add(v2)
        elif v2 == 6:
            connected_to_6.add(v1)
    
    # Color vertices connected to 6 with color 2
    for vertex in connected_to_6:
        coloring[vertex] = 2
        uncolored.remove(vertex)
    
    # Color remaining vertices
    for vertex in uncolored:
        # Try each color
        for color in colors:
            if is_safe(vertex, color, coloring, edges):
                coloring[vertex] = color
                break
    
    # Verify final solution
    for v1, v2 in edges:
        if coloring[v1] == coloring[v2]:
            return None
    
    return coloring

# Problem setup
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (0, 6), (2, 6), (2, 8), (3, 4), (3, 6), (5, 6), (5, 7)]
colors = [1, 2, 3]

# Find solution
solution = solve_coloring(vertices, edges, colors)

# Format and print solution
if solution:
    result = {str(k): v for k, v in solution.items()}
    print(result)