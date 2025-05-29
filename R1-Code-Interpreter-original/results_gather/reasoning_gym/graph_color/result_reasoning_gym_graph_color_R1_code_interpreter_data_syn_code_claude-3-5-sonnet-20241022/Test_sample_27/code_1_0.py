def is_valid(vertex, color, colored, edges):
    for v1, v2 in edges:
        if (v1 == vertex and v2 in colored and colored[v2] == color) or \
           (v2 == vertex and v1 in colored and colored[v1] == color):
            return False
    return True

def color_graph_backtrack(vertices, edges, colors):
    def backtrack(remaining_vertices):
        if not remaining_vertices:  # All vertices colored
            return colored.copy()
        
        vertex = remaining_vertices[0]
        remaining = remaining_vertices[1:]
        
        for color in colors:
            if is_valid(vertex, color, colored, edges):
                colored[vertex] = color
                result = backtrack(remaining)
                if result is not None:
                    return result
                del colored[vertex]
        return None

    # Start with vertex 6 (most connected)
    vertices_ordered = [6]  # Start with vertex 6
    # Add vertices connected to 6
    connected_to_6 = set()
    for v1, v2 in edges:
        if v1 == 6:
            connected_to_6.add(v2)
        elif v2 == 6:
            connected_to_6.add(v1)
    vertices_ordered.extend(sorted(connected_to_6))
    # Add remaining vertices
    remaining = set(vertices) - {6} - connected_to_6
    vertices_ordered.extend(sorted(remaining))

    colored = {}
    return backtrack(vertices_ordered)

# Problem setup
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 1), (0, 6), (2, 6), (2, 8), (3, 4), (3, 6), (5, 6), (5, 7)]
colors = [1, 2, 3]

# Find solution
solution = color_graph_backtrack(vertices, edges, colors)

# Format and print solution
if solution:
    result = {str(k): v for k, v in sorted(solution.items())}
    print(result)