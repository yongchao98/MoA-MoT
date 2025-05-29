def is_valid_coloring(graph, colors, vertex_colors):
    # Check if the current coloring is valid
    for v1, v2 in graph:
        if v1 in vertex_colors and v2 in vertex_colors:
            if vertex_colors[v1] == vertex_colors[v2]:
                return False
    return True

def find_coloring():
    # Graph definition
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(1, 9), (2, 3), (3, 7), (5, 9)]
    colors = [1, 2, 3]
    
    # Initialize all vertices with color 1
    vertex_colors = {v: 1 for v in vertices}
    
    # For each edge, ensure connected vertices have different colors
    for v1, v2 in edges:
        # If colors are same, change color of second vertex
        if vertex_colors[v1] == vertex_colors[v2]:
            # Try each color for v2 until we find a valid one
            for color in colors:
                vertex_colors[v2] = color
                if is_valid_coloring(edges, colors, vertex_colors):
                    break
    
    # Convert to JSON string format
    result = str(vertex_colors).replace("'", '"')
    print(f"<<<{result}>>>")

find_coloring()