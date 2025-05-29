import json

def is_safe(graph, vertex, color, coloring):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring():
    # Initialize graph
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 8), (2, 9), (3, 5), (6, 7), (7, 8)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    
    # Initialize coloring with first valid color for each vertex
    coloring = {}
    for vertex in vertices:
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                break
    
    # Convert to required format and print
    result = {str(k): v for k, v in coloring.items()}
    print(f"<<<{json.dumps(result)}>>>")

# Run the solution
graph_coloring()