import json

def is_safe(graph, vertex, colors, color, V):
    # Check if it's safe to color vertex with given color
    for i in range(V):
        if graph[vertex][i] == 1 and colors[i] == color:
            return False
    return True

def graph_coloring_util(graph, m, colors, vertex, V):
    # Base case: if all vertices are colored
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
V = 10  # number of vertices
m = 3   # number of colors
edges = [(0, 1), (3, 8), (4, 7), (4, 8), (4, 9), (8, 9)]

# Create adjacency matrix
graph = [[0 for _ in range(V)] for _ in range(V)]
for edge in edges:
    graph[edge[0]][edge[1]] = 1
    graph[edge[1]][edge[0]] = 1

# Initialize colors array
colors = [0] * V

# Find solution
if graph_coloring_util(graph, m, colors, 0, V):
    # Convert solution to required format
    solution = {str(i): colors[i] for i in range(V)}
    print(f"<<<{json.dumps(solution)}>>>")
else:
    print("No solution exists")