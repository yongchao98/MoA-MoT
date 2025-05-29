def is_valid_color(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

# Create adjacency list representation
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (0, 6), (1, 2), (3, 9)]
colors = [1, 2, 3]

# Initialize graph as adjacency list
graph = {v: [] for v in vertices}
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Initialize coloring
coloring = {}

# Color each vertex with minimum possible color
for vertex in vertices:
    for color in colors:
        if is_valid_color(graph, vertex, color, coloring):
            coloring[vertex] = color
            break

# Convert to string format and print
result = "{"
for v in vertices:
    result += f"\"{v}\": {coloring[v]}"
    if v != vertices[-1]:
        result += ", "
result += "}"
print(result)