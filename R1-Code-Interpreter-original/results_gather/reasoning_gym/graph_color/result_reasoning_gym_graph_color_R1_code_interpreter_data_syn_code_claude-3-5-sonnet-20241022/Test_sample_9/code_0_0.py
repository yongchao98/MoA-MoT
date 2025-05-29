import json

def is_safe(graph, vertex, colors, color, V):
    # Check if the color can be assigned to vertex
    for i in range(V):
        if graph[vertex][i] == 1 and colors[i] == color:
            return False
    return True

def graph_coloring_util(graph, m, colors, vertex, V):
    # Base case: If all vertices are colored
    if vertex == V:
        return True

    # Try different colors for vertex
    for color in range(1, m + 1):
        if is_safe(graph, vertex, colors, color, V):
            colors[vertex] = color
            if graph_coloring_util(graph, m, colors, vertex + 1, V):
                return True
            colors[vertex] = 0

    return False

# Initialize the problem
V = 10  # Number of vertices
m = 3   # Number of colors
edges = [(1, 2), (1, 5), (3, 5), (3, 9), (4, 5), (5, 7), (6, 8), (7, 9), (8, 9)]

# Create adjacency matrix
graph = [[0 for _ in range(V)] for _ in range(V)]
for edge in edges:
    graph[edge[0]][edge[1]] = 1
    graph[edge[1]][edge[0]] = 1

# Initialize colors array
colors = [0] * V

# Solve the coloring problem
if graph_coloring_util(graph, m, colors, 0, V):
    # Convert to required format
    result = {str(i): colors[i] for i in range(V)}
    print(f"<<<{json.dumps(result)}>>>")
else:
    print("No solution exists")