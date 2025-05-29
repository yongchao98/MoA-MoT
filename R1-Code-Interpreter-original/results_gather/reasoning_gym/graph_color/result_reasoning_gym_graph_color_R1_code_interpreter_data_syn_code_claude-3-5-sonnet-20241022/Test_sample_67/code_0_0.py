def is_valid_color(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, vertices, available_colors):
    colors = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in available_colors:
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                colors.pop(vertex, None)
        return False
    
    # Create adjacency list representation
    adj_list = {v: [] for v in vertices}
    edges = [(0, 5), (0, 9), (2, 7), (3, 4), (3, 6), (4, 6), (5, 6)]
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Try to color the graph
    if backtrack(0):
        # Convert all values to integers and ensure all vertices are included
        result = {str(v): colors.get(v, 1) for v in vertices}
        import json
        print(json.dumps(result))
    else:
        print("No solution found")

# Execute with given inputs
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
colors = [1, 2, 3]
graph_coloring(adj_list, vertices, colors)