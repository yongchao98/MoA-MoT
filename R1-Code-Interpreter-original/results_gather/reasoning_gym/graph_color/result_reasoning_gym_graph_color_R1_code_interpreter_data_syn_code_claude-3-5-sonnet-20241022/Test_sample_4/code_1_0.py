import json

def is_safe(adj_list, vertex, colors, color):
    # Check if any adjacent vertex has the same color
    for neighbor in adj_list[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def find_coloring():
    vertices = list(range(10))
    edges = [(0, 1), (0, 2), (0, 6), (1, 5), (3, 8), (6, 9)]
    num_colors = 3
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    colors = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
            
        current_vertex = vertices[vertex_index]
        for color in range(1, num_colors + 1):
            if is_safe(adj_list, current_vertex, colors, color):
                colors[current_vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                colors[current_vertex] = None
        return False
    
    if backtrack(0):
        # Convert the solution to the required format
        result = {str(k): v for k, v in colors.items()}
        print(json.dumps(result))
    else:
        print("No solution exists")

find_coloring()