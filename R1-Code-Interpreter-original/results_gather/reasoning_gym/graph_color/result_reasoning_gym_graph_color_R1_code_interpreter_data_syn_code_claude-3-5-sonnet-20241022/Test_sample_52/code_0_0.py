def is_valid_color(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(vertices, edges, colors):
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    
    # Initialize coloring dictionary
    coloring = {}
    
    # Try coloring each vertex
    for vertex in vertices:
        for color in colors:
            if is_valid_color(graph, vertex, color, coloring):
                coloring[vertex] = color
                break
    
    # Convert to string format for JSON
    result = {str(k): v for k, v in coloring.items()}
    print(result)

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 2), (0, 3), (2, 5), (3, 5), (4, 9)]
colors = [1, 2, 3]

graph_coloring(vertices, edges, colors)