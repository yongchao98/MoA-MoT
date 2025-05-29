def is_valid_coloring(graph, colors):
    # Check if adjacent vertices have different colors
    for v1, v2 in graph:
        if colors[v1] == colors[v2]:
            return False
    return True

def find_coloring():
    # Graph edges
    edges = [(0, 9), (2, 4), (2, 9), (5, 7)]
    vertices = list(range(10))
    possible_colors = [1, 2, 3]
    
    # Initialize all vertices with color 1
    coloring = {v: 1 for v in vertices}
    
    # For each edge, ensure connected vertices have different colors
    for v1, v2 in edges:
        if coloring[v1] == coloring[v2]:
            # Try to assign a different color to v2
            for color in possible_colors:
                if color != coloring[v1]:
                    coloring[v2] = color
                    break
    
    # Convert to JSON string format
    result = "{"
    for v in vertices:
        result += f'"{v}": {coloring[v]}'
        if v < vertices[-1]:
            result += ", "
    result += "}"
    
    print(f"<<<{result}>>>")

find_coloring()