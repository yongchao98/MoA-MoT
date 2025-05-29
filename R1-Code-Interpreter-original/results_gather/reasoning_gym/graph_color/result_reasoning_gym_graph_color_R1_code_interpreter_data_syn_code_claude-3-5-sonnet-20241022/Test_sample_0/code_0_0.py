def is_safe(graph, vertex, colors, color, V):
    # Check if the color can be assigned to vertex
    for i in range(V):
        if graph[vertex][i] == 1 and colors[i] == color:
            return False
    return True

def graph_coloring_util(graph, m, colors, vertex, V):
    # Base case: all vertices are colored
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

# Initialize the graph
V = 10  # vertices
graph = [[0] * V for _ in range(V)]
edges = [(0, 5), (0, 6), (0, 9), (2, 9), (5, 9), (6, 7), (7, 8)]

# Create adjacency matrix
for edge in edges:
    graph[edge[0]][edge[1]] = 1
    graph[edge[1]][edge[0]] = 1

# Initialize colors array
colors = [0] * V

# Number of colors
m = 3

# Solve the graph coloring problem
if graph_coloring_util(graph, m, colors, 0, V):
    # Convert to required JSON format
    result = {str(i): colors[i] for i in range(V)}
    print(result)
else:
    print("No solution exists")