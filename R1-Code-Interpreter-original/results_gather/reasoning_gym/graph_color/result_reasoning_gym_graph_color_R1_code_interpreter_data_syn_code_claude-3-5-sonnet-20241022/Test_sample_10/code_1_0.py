def is_valid_color(adj_list, colors, vertex, color):
    for neighbor in adj_list[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def find_coloring():
    # Initialize adjacency list
    vertices = list(range(10))
    edges = [(0, 1), (0, 7), (0, 9), (1, 4), (2, 4), (3, 5), (3, 6), (6, 8), (7, 9)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for u, v in edges:
        adj_list[u].append(v)
        adj_list[v].append(u)
    
    # Initialize colors dictionary
    vertex_colors = {}
    
    def backtrack(vertex):
        if vertex == len(vertices):
            return True
            
        for color in colors:
            if is_valid_color(adj_list, vertex_colors, vertex, color):
                vertex_colors[vertex] = color
                if backtrack(vertex + 1):
                    return True
                vertex_colors.pop(vertex)
        return False
    
    # Find solution
    backtrack(0)
    
    # Convert to string keys for JSON format
    result = {str(k): v for k, v in vertex_colors.items()}
    print(f"<<<{result}>>>")

find_coloring()