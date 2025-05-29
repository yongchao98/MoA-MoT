def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(vertices, edges, colors):
    # Initialize graph as adjacency list
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    
    # Initialize coloring dictionary
    coloring = {}
    
    # Try to color each vertex
    for vertex in vertices:
        # Try each color
        colored = False
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                colored = True
                break
        # If no color works, try the first available color
        if not colored:
            coloring[vertex] = colors[0]
    
    # Convert result to JSON string
    result = {str(k): v for k, v in coloring.items()}
    print("<<<" + str(result) + ">>>")

# Run the solution
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 9), (2, 3), (2, 6), (3, 8), (5, 7)]
colors = [1, 2, 3]
graph_coloring(vertices, edges, colors)