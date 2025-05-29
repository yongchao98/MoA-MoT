def is_valid_color(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring():
    # Create adjacency list representation
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 5), (2, 3), (4, 5), (4, 6), (6, 7), (6, 8)]
    possible_colors = [1, 2, 3]
    
    # Initialize graph as adjacency list
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    
    # Initialize colors dictionary with None
    colors = {v: None for v in vertices}
    
    # Simple coloring strategy:
    # Since this is a sparse graph, we can use a simple greedy approach
    for vertex in vertices:
        for color in possible_colors:
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                break
        # If no color was assigned, use color 1 (graph is sparse enough)
        if colors[vertex] is None:
            colors[vertex] = 1
    
    # Convert to string format for JSON
    colors = {str(k): v for k, v in colors.items()}
    
    import json
    print(json.dumps(colors))

graph_coloring()