def is_valid_color(graph, colors, vertex, color):
    for v in range(len(graph)):
        if graph[vertex][v] == 1 and colors[v] == color:
            return False
    return True

def graph_coloring(graph, m, v, colors, V):
    if v == V:
        return True
    
    for c in range(1, m + 1):
        if is_valid_color(graph, colors, v, c):
            colors[v] = c
            if graph_coloring(graph, m, v + 1, colors, V):
                return True
            colors[v] = 0
    return False

# Create adjacency matrix
V = 10  # number of vertices
graph = [[0] * V for _ in range(V)]
edges = [(0, 7), (0, 9), (1, 4), (2, 8), (4, 5), (7, 8), (8, 9)]

# Fill adjacency matrix
for edge in edges:
    graph[edge[0]][edge[1]] = 1
    graph[edge[1]][edge[0]] = 1

# Initialize colors array
colors = [0] * V
m = 3  # number of colors

# Find solution
if graph_coloring(graph, m, 0, colors, V):
    # Convert to required format
    result = {str(i): colors[i] for i in range(V)}
    print(result)
else:
    print("No solution exists")