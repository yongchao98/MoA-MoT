import json

def is_safe(graph, vertex, colors, color):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, vertices):
    colors = {}
    for vertex in vertices:
        colors[vertex] = None
    
    def backtrack(vertex):
        if vertex not in vertices:
            return True
        
        for color in range(1, num_colors + 1):
            if is_safe(graph, vertex, colors, color):
                colors[vertex] = color
                next_vertex = vertices[vertices.index(vertex) + 1] if vertices.index(vertex) < len(vertices) - 1 else None
                if next_vertex is None or backtrack(next_vertex):
                    return True
                colors[vertex] = None
        return False

    # Create adjacency list representation
    adj_list = {v: [] for v in vertices}
    edges = [(0, 8), (1, 9), (3, 6), (4, 8), (5, 7), (5, 9)]
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Start with vertex 0
    if backtrack(vertices[0]):
        # Convert all values to integers and create the result
        result = {str(k): v for k, v in colors.items()}
        print(f"<<<{json.dumps(result)}>>>")
    else:
        print("No solution exists")

# Set up the problem
vertices = list(range(10))
num_colors = 3

# Solve the problem
graph_coloring(adj_list, num_colors, vertices)