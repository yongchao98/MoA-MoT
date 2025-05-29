def is_valid_color(edges, colors, vertex, color):
    for v1, v2 in edges:
        if v1 == vertex and colors[v2] == color:
            return False
        if v2 == vertex and colors[v1] == color:
            return False
    return True

def graph_coloring():
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 5), (0, 7), (2, 3), (2, 7), (3, 5), (4, 7)]
    possible_colors = [1, 2, 3]
    colors = [-1] * len(vertices)
    
    def backtrack(vertex):
        if vertex == len(vertices):
            return True
            
        for color in possible_colors:
            if is_valid_color(edges, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex + 1):
                    return True
                colors[vertex] = -1
        return False
    
    if backtrack(0):
        result = {str(v): colors[v] for v in vertices}
        print(str(result).replace("'", '"'))
    else:
        print("No solution exists")

graph_coloring()