import json

def is_valid_color(graph, colors, vertex, color):
    # Check if any neighbor has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_vertices, available_colors):
    colors = {}
    
    def backtrack(vertex):
        if vertex == num_vertices:
            return True
        
        # Try each color for current vertex
        for color in available_colors:
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex + 1):
                    return True
                colors.pop(vertex, None)
        
        return False
    
    # Create adjacency list representation
    adj_list = {i: [] for i in range(num_vertices)}
    edges = [(0, 2), (0, 5), (1, 2), (1, 9), (2, 6), (5, 9)]
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Find coloring
    backtrack(0)
    
    # Convert all values to integers and ensure all vertices are colored
    result = {str(i): colors.get(i, 1) for i in range(num_vertices)}
    print(json.dumps(result))

# Run the coloring algorithm
graph_coloring({i: [] for i in range(10)}, 10, [1, 2, 3])