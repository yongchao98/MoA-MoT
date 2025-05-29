import json

def is_safe(graph, vertex, color, coloring):
    # Check if it's safe to color vertex with given color
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
    
    coloring = {}
    
    # Start with first vertex
    for v in vertices:
        if v not in coloring:
            # Try each color
            for color in colors:
                if is_safe(graph, v, color, coloring):
                    coloring[v] = color
                    break
    
    # Convert all values to strings for JSON
    result = {str(k): v for k, v in coloring.items()}
    print(f"<<<{json.dumps(result)}>>>")

# Define the problem parameters
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 5), (1, 8), (3, 4), (3, 5), (5, 6)]
colors = [1, 2, 3]

# Solve the problem
graph_coloring(vertices, edges, colors)