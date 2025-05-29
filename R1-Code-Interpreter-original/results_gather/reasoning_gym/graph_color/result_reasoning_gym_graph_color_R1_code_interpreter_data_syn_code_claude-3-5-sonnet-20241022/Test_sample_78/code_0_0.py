import json

def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

# Create graph representation (adjacency list)
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (1, 6), (4, 7), (7, 8)]

# Create adjacency list
graph = {v: [] for v in vertices}
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Initialize colors dictionary
colors = {}

# Color vertices (simple greedy approach)
for vertex in vertices:
    # Try each color
    for color in [1, 2, 3]:
        if is_valid_coloring(graph, colors, vertex, color):
            colors[vertex] = color
            break

# Convert to string format required
result = {str(k): v for k, v in colors.items()}
print(f"<<<{json.dumps(result)}>>>")