def graph_coloring(vertices, edges, colors):
    # Initialize the color map with None
    color_map = {v: None for v in vertices}
    
    # Function to check if the current color assignment is safe for vertex v
    def is_safe(v, c):
        for edge in edges:
            if v in edge:
                # Get the adjacent vertex
                adjacent = edge[0] if edge[1] == v else edge[1]
                if color_map[adjacent] == c:
                    return False
        return True
    
    # Recursive function to solve the coloring problem
    def solve(v_index):
        if v_index == len(vertices):
            return True
        
        vertex = vertices[v_index]
        for color in colors:
            if is_safe(vertex, color):
                color_map[vertex] = color
                if solve(v_index + 1):
                    return True
                # Backtrack
                color_map[vertex] = None
        
        return False
    
    # Start solving from the first vertex
    solve(0)
    return color_map

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 8), (1, 4), (1, 7), (1, 8), (8, 9)]
colors = [1, 2, 3]

# Get the coloring solution
solution = graph_coloring(vertices, edges, colors)
print(solution)