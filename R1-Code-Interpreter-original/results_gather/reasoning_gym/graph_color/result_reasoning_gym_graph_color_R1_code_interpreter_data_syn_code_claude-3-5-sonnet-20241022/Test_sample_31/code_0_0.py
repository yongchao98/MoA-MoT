def is_valid_color(graph, colors, vertex, color):
    for neighbor in range(len(graph)):
        if graph[vertex][neighbor] == 1 and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, m, colors, vertex, V):
    if vertex == V:
        return True
    
    for color in range(1, m + 1):
        if is_valid_color(graph, colors, vertex, color):
            colors[vertex] = color
            if graph_coloring(graph, m, colors, vertex + 1, V):
                return True
            colors[vertex] = 0
    return False

# Create adjacency matrix from edges
V = 10  # number of vertices
graph = [[0] * V for _ in range(V)]
edges = [(0, 6), (1, 4), (2, 3), (2, 4), (3, 8), (4, 8), (5, 7), (6, 7), (6, 8), (6, 9)]
for edge in edges:
    graph[edge[0]][edge[1]] = 1
    graph[edge[1]][edge[0]] = 1

# Initialize colors array
colors = [0] * V

# Try coloring with 3 colors
m = 3
if graph_coloring(graph, m, colors, 0, V):
    # Convert result to required JSON format
    result = {str(i): colors[i] for i in range(V)}
    print(f"<<<{result}>>>")
else:
    print("No solution exists")