def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring():
    # Create adjacency list representation
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 8), (1, 4), (1, 7), (1, 8), (8, 9)]
    colors = [1, 2, 3]
    
    # Build adjacency list
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    coloring = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                coloring.pop(vertex, None)
        return False

    # Find solution using backtracking
    if backtrack(0):
        # Add isolated vertices with color 1
        for v in vertices:
            if v not in coloring:
                coloring[v] = 1
        # Convert to JSON-formatted string
        result = "<<<" + str(coloring).replace("'", '"') + ">>>"
        print(result)
    else:
        print("No solution exists")

# Run the solution
graph_coloring()