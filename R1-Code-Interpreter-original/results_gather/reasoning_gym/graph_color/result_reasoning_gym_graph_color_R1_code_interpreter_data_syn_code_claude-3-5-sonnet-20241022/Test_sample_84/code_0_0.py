def create_adjacency_list(vertices, edges):
    adj_list = {v: [] for v in vertices}
    for u, v in edges:
        adj_list[u].append(v)
        adj_list[v].append(u)
    return adj_list

def is_valid_color(vertex, color, coloring, adj_list):
    for neighbor in adj_list[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def color_graph(vertices, edges, colors):
    adj_list = create_adjacency_list(vertices, edges)
    coloring = {}
    
    # For each vertex
    for vertex in vertices:
        # Try each color
        for color in colors:
            if is_valid_color(vertex, color, coloring, adj_list):
                coloring[vertex] = color
                break
    
    # Convert to string format for JSON
    result = {str(k): v for k, v in coloring.items()}
    print(result)

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (1, 7), (2, 8), (5, 9)]
colors = [1, 2, 3]

color_graph(vertices, edges, colors)