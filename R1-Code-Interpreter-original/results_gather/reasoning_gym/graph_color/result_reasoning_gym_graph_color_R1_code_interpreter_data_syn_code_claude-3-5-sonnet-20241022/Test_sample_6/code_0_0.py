import json

def is_valid_coloring(graph, colors, vertex_colors):
    # Check if any adjacent vertices have the same color
    for v1, v2 in graph:
        if v1 in vertex_colors and v2 in vertex_colors:
            if vertex_colors[v1] == vertex_colors[v2]:
                return False
    return True

def color_graph():
    # Graph definition
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 3), (2, 5), (3, 5), (4, 5)]
    colors = [1, 2, 3]
    
    # Initialize all vertices with color 1
    vertex_colors = {v: 1 for v in vertices}
    
    # For vertices that are part of edges, assign minimum valid color
    connected_vertices = set()
    for v1, v2 in edges:
        connected_vertices.add(v1)
        connected_vertices.add(v2)
    
    # Color connected vertices
    for vertex in connected_vertices:
        for color in colors:
            vertex_colors[vertex] = color
            if is_valid_coloring(edges, colors, vertex_colors):
                break
    
    # Convert to JSON string format
    result = json.dumps(vertex_colors)
    print(result)

color_graph()