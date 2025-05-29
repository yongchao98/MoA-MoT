import json

def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

# Create graph representation as adjacency list
vertices = list(range(10))
edges = [(2, 4), (3, 7), (5, 8), (7, 9)]
possible_colors = [1, 2, 3]

# Convert edges to adjacency list
graph = {v: [] for v in vertices}
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Create coloring
colors = {}
for vertex in vertices:
    # Try each color
    for color in possible_colors:
        if is_valid_coloring(graph, colors, vertex, color):
            colors[vertex] = color
            break

# Convert to string format and print
result = {str(k): v for k, v in colors.items()}
print(json.dumps(result))