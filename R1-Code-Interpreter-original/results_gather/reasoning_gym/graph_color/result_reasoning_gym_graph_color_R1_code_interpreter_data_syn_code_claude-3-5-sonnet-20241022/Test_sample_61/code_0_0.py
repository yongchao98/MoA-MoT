def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, available_colors):
    n = len(graph)
    colors = {}
    
    def backtrack(vertex):
        if vertex == n:
            return True
        
        # Skip vertices that have no edges (can be any color)
        if vertex not in graph:
            colors[vertex] = available_colors[0]
            return backtrack(vertex + 1)
        
        # Try each color for the current vertex
        for color in available_colors:
            if is_valid_coloring(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex + 1):
                    return True
        
        return False
    
    # Create adjacency list representation
    adj_list = {}
    edges = [(0, 4), (1, 7), (3, 5), (3, 6), (4, 5), (6, 8)]
    for v1, v2 in edges:
        if v1 not in adj_list:
            adj_list[v1] = []
        if v2 not in adj_list:
            adj_list[v2] = []
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize colors for isolated vertices
    for vertex in range(n):
        if vertex not in adj_list:
            colors[vertex] = available_colors[0]
    
    backtrack(0)
    
    # Convert all values to strings for JSON format
    result = {str(k): v for k, v in colors.items()}
    print(result)

# Run the coloring algorithm
vertices = list(range(10))
colors = [1, 2, 3]
graph_coloring(vertices, colors)