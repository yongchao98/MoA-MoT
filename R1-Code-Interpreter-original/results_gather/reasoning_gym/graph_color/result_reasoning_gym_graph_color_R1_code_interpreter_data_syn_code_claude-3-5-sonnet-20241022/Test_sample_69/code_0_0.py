import json

def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

# Create graph representation as adjacency list
vertices = list(range(10))
edges = [(0, 4), (2, 6), (4, 5), (5, 6)]

# Convert to adjacency list
graph = {v: [] for v in vertices}
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Initialize colors dictionary
colors = {}

# For each vertex, assign the first valid color
for vertex in vertices:
    for color in [1, 2, 3]:
        if is_valid_coloring(graph, colors, vertex, color):
            colors[vertex] = color
            break

# Convert to string format and print
result = {str(k): v for k, v in colors.items()}
print(f"<<<{json.dumps(result)}>>>")