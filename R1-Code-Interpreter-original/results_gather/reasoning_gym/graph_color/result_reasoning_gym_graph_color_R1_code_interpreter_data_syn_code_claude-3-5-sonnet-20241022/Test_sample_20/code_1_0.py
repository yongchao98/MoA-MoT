import json

def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(vertices, edges, available_colors):
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    
    colors = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
            
        vertex = vertices[vertex_index]
        # Sort colors to try minimal colors first
        for color in sorted(available_colors):
            if is_valid_coloring(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                colors.pop(vertex, None)
        return False
    
    # Start with vertices that have most connections
    vertex_degrees = [(v, len(graph[v])) for v in vertices]
    vertices = [v for v, d in sorted(vertex_degrees, key=lambda x: x[1], reverse=True)]
    
    backtrack(0)
    
    # Convert all keys to strings for JSON format
    result = {str(k): v for k, v in colors.items()}
    print(json.dumps(result))

# Problem inputs
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (1, 5), (1, 7), (1, 8), (4, 8), (6, 7), (6, 8)]
available_colors = [1, 2, 3]

graph_coloring(vertices, edges, available_colors)