def is_valid_color(graph, colors, vertex, color):
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
    
    # Sort vertices by degree (number of connections)
    vertices_by_degree = sorted(vertices, key=lambda v: len(graph[v]), reverse=True)
    
    colors = {}
    
    def color_vertex(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices_by_degree[vertex_index]
        for color in available_colors:
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                if color_vertex(vertex_index + 1):
                    return True
                colors.pop(vertex, None)
        return False
    
    if color_vertex(0):
        # Convert to string keys for JSON compatibility
        result = {str(k): v for k, v in colors.items()}
        print(f"<<<%s>>>" % str(result).replace("'", '"'))
        return True
    else:
        print("No solution found")
        return False

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 7), (2, 4), (2, 8), (2, 9), (4, 6), (5, 8), (5, 9), (7, 9)]
available_colors = [1, 2, 3]

graph_coloring(vertices, edges, available_colors)