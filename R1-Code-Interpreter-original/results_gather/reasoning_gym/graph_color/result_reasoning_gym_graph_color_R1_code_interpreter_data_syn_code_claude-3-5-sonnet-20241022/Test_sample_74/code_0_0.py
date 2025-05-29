import json

# Create adjacency list representation
def create_graph():
    graph = {i: [] for i in range(10)}
    edges = [(0, 1), (0, 4), (2, 4), (2, 9)]
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    return graph

def get_valid_color(vertex, colored, graph):
    neighbor_colors = set(colored.get(neighbor) for neighbor in graph[vertex] if neighbor in colored)
    for color in [1, 2, 3]:
        if color not in neighbor_colors:
            return color
    return None

# Color the graph
def color_graph():
    graph = create_graph()
    colored = {}
    
    # Color each vertex
    for vertex in range(10):
        color = get_valid_color(vertex, colored, graph)
        colored[vertex] = color if color else 1  # If no constraints, use color 1
    
    # Convert to required format
    result = {str(k): v for k, v in colored.items()}
    print(json.dumps(result))

color_graph()