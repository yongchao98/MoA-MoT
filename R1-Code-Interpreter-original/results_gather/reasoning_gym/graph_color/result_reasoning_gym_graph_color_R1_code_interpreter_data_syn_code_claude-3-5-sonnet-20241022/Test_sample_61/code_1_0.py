def is_valid_coloring(adj_list, colors, vertex, color):
    if vertex not in adj_list:
        return True
    for neighbor in adj_list[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def find_coloring():
    # Graph properties
    vertices = list(range(10))
    edges = [(0, 4), (1, 7), (3, 5), (3, 6), (4, 5), (6, 8)]
    available_colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    colors = {}
    
    def backtrack(vertex):
        if vertex == len(vertices):
            return True
            
        for color in available_colors:
            if is_valid_coloring(adj_list, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex + 1):
                    return True
                colors.pop(vertex, None)
        return False
    
    # Find a valid coloring
    backtrack(0)
    
    # Convert to required format
    result = {str(k): v for k, v in colors.items()}
    print(result)

find_coloring()