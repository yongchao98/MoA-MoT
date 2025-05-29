import json

def is_safe(graph, vertex, colors, color):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, colors, vertex):
    # If all vertices are colored, return True
    if vertex == len(graph):
        return True
    
    # Try each color for current vertex
    for color in range(1, num_colors + 1):
        if is_safe(graph, vertex, colors, color):
            colors[vertex] = color
            if graph_coloring(graph, num_colors, colors, vertex + 1):
                return True
            colors[vertex] = 0
    
    return False

# Create adjacency list representation
V = 10
graph = {i: [] for i in range(V)}
edges = [(0, 4), (0, 6), (3, 4), (4, 9), (5, 7)]

# Build the graph
for u, v in edges:
    graph[u].append(v)
    graph[v].append(u)

# Initialize colors array
colors = [0] * V

# Find solution
graph_coloring(graph, 3, colors, 0)

# Convert to required format
result = {str(i): colors[i] for i in range(V)}

# Print result in required format
print(f"<<<{json.dumps(result)}>>>")