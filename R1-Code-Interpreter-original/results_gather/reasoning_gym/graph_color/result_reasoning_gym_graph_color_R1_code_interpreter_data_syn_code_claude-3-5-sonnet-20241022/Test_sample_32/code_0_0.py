import json

def is_safe(graph, vertex, color, coloring):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
    coloring = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
        return False
    
    # Create adjacency list representation
    graph_adj = {v: [] for v in vertices}
    for edge in [(1, 5), (1, 8), (3, 4), (3, 5), (5, 6)]:
        graph_adj[edge[0]].append(edge[1])
        graph_adj[edge[1]].append(edge[0])
    
    # Try to color the graph
    backtrack(0)
    
    # Convert all values to strings for JSON
    result = {str(k): v for k, v in coloring.items()}
    print(f"<<<{json.dumps(result)}>>>")

# Define the problem parameters
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
colors = [1, 2, 3]

# Solve the problem
graph_coloring({v: [] for v in vertices}, colors, vertices)