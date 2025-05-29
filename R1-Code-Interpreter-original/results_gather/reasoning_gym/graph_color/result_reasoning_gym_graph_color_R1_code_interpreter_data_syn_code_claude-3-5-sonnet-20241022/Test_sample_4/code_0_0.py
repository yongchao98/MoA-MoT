import json

def is_safe(graph, vertex, colors, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, vertices):
    colors = {}
    for vertex in vertices:
        colors[vertex] = None
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        current_vertex = vertices[vertex_index]
        for color in range(1, num_colors + 1):
            if is_safe(graph, current_vertex, colors, color):
                colors[current_vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                colors[current_vertex] = None
        return False

    # Create adjacency list representation
    adj_list = {v: [] for v in vertices}
    edges = [(0, 1), (0, 2), (0, 6), (1, 5), (3, 8), (6, 9)]
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)

    # Try to color the graph
    if backtrack(0):
        # Convert all values to integers and create the result dictionary
        result = {str(k): v for k, v in colors.items()}
        print(json.dumps(result))
    else:
        print("No solution exists")

# Execute with given parameters
vertices = list(range(10))
num_colors = 3
graph_coloring(vertices, num_colors, vertices)