import json

def is_valid_coloring(edges, vertex_colors):
    # Check if any adjacent vertices have the same color
    for v1, v2 in edges:
        if vertex_colors[v1] == vertex_colors[v2]:
            return False
    return True

def get_valid_color(vertex, edges, vertex_colors, colors):
    # Get neighbors of the vertex
    neighbors = set()
    for v1, v2 in edges:
        if v1 == vertex:
            neighbors.add(v2)
        elif v2 == vertex:
            neighbors.add(v1)
    
    # Find a color that's not used by any neighbor
    neighbor_colors = {vertex_colors[n] for n in neighbors if n in vertex_colors}
    for color in colors:
        if color not in neighbor_colors:
            return color
    return None

def color_graph():
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 3), (2, 5), (3, 5), (4, 5)]
    colors = [1, 2, 3]
    
    # Initialize all vertices with color 1
    vertex_colors = {v: 1 for v in vertices}
    
    # Get all vertices that are part of edges
    connected_vertices = set()
    for v1, v2 in edges:
        connected_vertices.add(v1)
        connected_vertices.add(v2)
    
    # Color connected vertices one by one
    for vertex in connected_vertices:
        color = get_valid_color(vertex, edges, vertex_colors, colors)
        if color:
            vertex_colors[vertex] = color
    
    # Verify the coloring is valid
    if is_valid_coloring(edges, vertex_colors):
        print(json.dumps(vertex_colors))

color_graph()