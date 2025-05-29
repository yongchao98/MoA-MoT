import json

def is_safe(graph, vertex, colors, color):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, vertex, colors):
    # Base case: if all vertices are colored
    if vertex == len(graph):
        return True
    
    # Try each color for current vertex
    for color in range(1, num_colors + 1):
        if is_safe(graph, vertex, colors, color):
            colors[vertex] = color
            if graph_coloring(graph, num_colors, vertex + 1, colors):
                return True
            colors[vertex] = 0
    
    return False

# Create adjacency list representation
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 5), (0, 9), (1, 9), (2, 7), (3, 6), (5, 7)]
num_colors = 3

# Initialize graph as adjacency list
graph = {v: [] for v in vertices}
for edge in edges:
    graph[edge[0]].append(edge[1])
    graph[edge[1]].append(edge[0])

# Initialize colors array
colors = [0] * len(vertices)

# Find solution
graph_coloring(graph, num_colors, 0, colors)

# Convert to required format
result = {str(v): colors[v] for v in vertices}
print(json.dumps(result))