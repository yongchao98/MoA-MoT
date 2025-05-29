def is_valid_color(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(vertices, edges, available_colors):
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    
    # Initialize colors dictionary
    colors = {}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in available_colors:
            if is_valid_color(graph, colors, vertex, color):
                colors[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
        return False
    
    # Start backtracking from first vertex
    backtrack(0)
    
    # Convert all values to strings for JSON compatibility
    result = {str(k): v for k, v in colors.items()}
    print(result)

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 3), (0, 8), (2, 5), (4, 6)]
colors = [1, 2, 3]

graph_coloring(vertices, edges, colors)