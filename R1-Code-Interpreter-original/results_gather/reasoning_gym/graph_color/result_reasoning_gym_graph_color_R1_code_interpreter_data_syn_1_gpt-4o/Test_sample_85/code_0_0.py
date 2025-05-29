def is_safe(vertex, color, graph, colors):
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, m, colors, vertex):
    if vertex == len(graph):
        return True

    for color in range(1, m + 1):
        if is_safe(vertex, color, graph, colors):
            colors[vertex] = color
            if graph_coloring(graph, m, colors, vertex + 1):
                return True
            colors[vertex] = 0

    return False

def solve_graph_coloring():
    vertices = 10
    edges = [(1, 5), (3, 8), (4, 5), (7, 8), (8, 9)]
    graph = {i: [] for i in range(vertices)}
    for u, v in edges:
        graph[u].append(v)
        graph[v].append(u)

    colors = [0] * vertices
    m = 3  # Number of colors

    if graph_coloring(graph, m, colors, 0):
        return {str(i): colors[i] for i in range(vertices)}
    else:
        return "No solution"

result = solve_graph_coloring()
print(result)