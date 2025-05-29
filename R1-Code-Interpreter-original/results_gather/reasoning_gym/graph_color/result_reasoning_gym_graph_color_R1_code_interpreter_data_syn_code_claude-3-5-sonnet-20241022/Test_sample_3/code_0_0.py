import json

def create_graph(vertices, edges):
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    return graph

def is_safe(graph, vertex, color, coloring):
    # Check if it's safe to color the vertex with given color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
    coloring = {}
    
    def backtrack(vertex_idx):
        if vertex_idx == len(vertices):
            return True
        
        vertex = vertices[vertex_idx]
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                if backtrack(vertex_idx + 1):
                    return True
                coloring.pop(vertex, None)
        return False
    
    # Start backtracking from first vertex
    backtrack(0)
    return coloring

# Problem inputs
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 9), (2, 5), (2, 7), (3, 6), (3, 8), (4, 5)]
colors = [1, 2, 3]

# Create graph and find coloring
graph = create_graph(vertices, edges)
coloring = graph_coloring(graph, colors, vertices)

# Convert to string format required
result = json.dumps({str(k): v for k, v in coloring.items()})
print(result)