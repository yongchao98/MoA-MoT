import json

def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

# Create graph representation (adjacency list)
graph = {i: [] for i in range(10)}
edges = [(2, 4), (3, 9), (4, 9), (6, 7)]
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Initialize colors dictionary
colors = {}

# Try coloring each vertex
available_colors = [1, 2, 3]
for vertex in range(10):
    # Try each color
    for color in available_colors:
        if is_valid_coloring(graph, colors, vertex, color):
            colors[vertex] = color
            break

# Convert all values to strings for JSON format
result = {str(k): v for k, v in colors.items()}

# Print the result in required format
print(f"<<<{json.dumps(result)}>>>")