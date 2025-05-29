import json

def is_safe(graph, vertex, colors, color):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, colors, vertex):
    # Base case: if all vertices are colored
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

# Create graph representation (adjacency list)
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (1, 2), (1, 3), (1, 9), (3, 5), (4, 8), (7, 9)]

# Create adjacency list
graph = {v: [] for v in vertices}
for edge in edges:
    graph[edge[0]].append(edge[1])
    graph[edge[1]].append(edge[0])

# Initialize colors array
colors = [0] * len(vertices)

# Find solution
graph_coloring(graph, 3, colors, 0)

# Create solution dictionary
solution = {str(v): colors[v] for v in vertices}

# Print solution in required format
print(f"<<<{json.dumps(solution)}>>>")