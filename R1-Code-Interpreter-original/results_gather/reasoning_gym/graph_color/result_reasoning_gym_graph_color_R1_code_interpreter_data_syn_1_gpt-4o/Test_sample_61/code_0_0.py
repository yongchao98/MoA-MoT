def is_valid_color(vertex, color, graph, colors):
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, colors, vertex):
    if vertex == len(graph):
        return True
    
    for color in range(1, num_colors + 1):
        if is_valid_color(vertex, color, graph, colors):
            colors[vertex] = color
            if graph_coloring(graph, num_colors, colors, vertex + 1):
                return True
            colors[vertex] = 0
    
    return False

def solve_graph_coloring():
    vertices = 10
    edges = [(0, 4), (1, 7), (3, 5), (3, 6), (4, 5), (6, 8)]
    num_colors = 3
    
    # Create adjacency list for the graph
    graph = {i: [] for i in range(vertices)}
    for (u, v) in edges:
        graph[u].append(v)
        graph[v].append(u)
    
    # Initialize colors for each vertex
    colors = [0] * vertices
    
    if graph_coloring(graph, num_colors, colors, 0):
        return {str(i): colors[i] for i in range(vertices)}
    else:
        return "No solution"

# Execute the function and print the result
result = solve_graph_coloring()
print(result)