import json

def build_graph(vertices, edges):
    # Create adjacency list
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    return graph

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

    coloring = {}
    if backtrack(0):
        return coloring
    return None

# Problem inputs
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (1, 2), (1, 6), (4, 6), (5, 6)]
colors = [1, 2, 3]

# Solve the problem
graph = build_graph(vertices, edges)
solution = graph_coloring(graph, colors, vertices)

# Convert solution to string format
result = {str(k): v for k, v in solution.items()}
print(json.dumps(result))