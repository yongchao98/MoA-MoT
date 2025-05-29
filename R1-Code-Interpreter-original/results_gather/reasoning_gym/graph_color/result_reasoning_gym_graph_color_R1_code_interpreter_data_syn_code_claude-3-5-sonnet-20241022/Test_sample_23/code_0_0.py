def is_valid_color(graph, colors, vertex, color):
    for neighbor in range(len(graph)):
        if graph[vertex][neighbor] == 1 and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, vertices):
    colors = [-1] * len(vertices)
    
    def backtrack(vertex):
        if vertex == len(vertices):
            return True
            
        for color in range(1, num_colors + 1):
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex + 1):
                    return True
                colors[vertex] = -1
        return False
    
    # Create adjacency matrix
    n = len(vertices)
    adj_matrix = [[0] * n for _ in range(n)]
    edges = [(0, 5), (0, 7), (2, 3), (2, 7), (3, 5), (4, 7)]
    for v1, v2 in edges:
        adj_matrix[v1][v2] = 1
        adj_matrix[v2][v1] = 1
    
    # Find solution
    if backtrack(0):
        # Convert result to JSON format
        result = {str(v): colors[v] for v in vertices}
        print(str(result).replace("'", '"'))
    else:
        print("No solution exists")

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
num_colors = 3

graph_coloring(graph=[[0]*10 for _ in range(10)], num_colors=num_colors, vertices=vertices)