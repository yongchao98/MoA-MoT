def is_valid_color(graph, colors, vertex, color):
    for neighbor in range(len(graph)):
        if graph[vertex][neighbor] == 1 and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, m, vertex, colors, V):
    if vertex == V:
        return True
    
    for color in range(1, m + 1):
        if is_valid_color(graph, colors, vertex, color):
            colors[vertex] = color
            if graph_coloring(graph, m, vertex + 1, colors, V):
                return True
            colors[vertex] = 0
    return False

# Create adjacency matrix from edges
V = 10  # number of vertices
edges = [(0, 1), (0, 8), (1, 4), (1, 6), (2, 7), (3, 7), (4, 9)]
graph = [[0] * V for _ in range(V)]
for u, v in edges:
    graph[u][v] = graph[v][u] = 1

# Initialize colors array
colors = [0] * V

# Try coloring with 3 colors
m = 3
if graph_coloring(graph, m, 0, colors, V):
    # Convert result to JSON format
    result = {str(i): colors[i] for i in range(V)}
    print(result)
else:
    print("No solution exists")