import json

def is_valid_color(adj_list, colors, vertex, color):
    # Check if any neighbor has the same color
    for neighbor in adj_list[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring():
    # Initialize vertices and edges
    vertices = list(range(10))
    edges = [(0, 2), (0, 5), (1, 2), (1, 9), (2, 6), (5, 9)]
    colors = [1, 2, 3]
    
    # Create adjacency list
    adj_list = {v: [] for v in vertices}
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # Initialize coloring with first color
    vertex_colors = {}
    
    def color_vertex(vertex):
        if vertex == len(vertices):
            return True
            
        for color in colors:
            if is_valid_color(adj_list, vertex_colors, vertex, color):
                vertex_colors[vertex] = color
                if color_vertex(vertex + 1):
                    return True
                vertex_colors.pop(vertex, None)
        return False
    
    # Start coloring from vertex 0
    color_vertex(0)
    
    # Fill in any uncolored vertices with color 1 (for isolated vertices)
    for v in vertices:
        if v not in vertex_colors:
            vertex_colors[v] = 1
    
    # Convert to required format
    result = {str(v): vertex_colors[v] for v in vertices}
    print(json.dumps(result))

graph_coloring()