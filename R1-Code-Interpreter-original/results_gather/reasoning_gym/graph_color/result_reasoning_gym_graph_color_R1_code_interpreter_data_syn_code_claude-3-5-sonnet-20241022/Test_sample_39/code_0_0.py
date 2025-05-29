def is_safe(graph, vertex, color, coloring):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
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

    # Create adjacency list representation
    adj_list = {v: [] for v in vertices}
    edges = [(0, 7), (1, 9), (3, 4), (3, 6), (5, 6), (5, 7), (5, 8)]
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)

    # Initialize coloring dictionary
    coloring = {}
    colors = [1, 2, 3]
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    # Find solution using backtracking
    backtrack(0)

    # Convert result to string format
    result = {str(k): v for k, v in coloring.items()}
    print(str(result).replace("'", '"'))

graph_coloring({}, [1, 2, 3], range(10))