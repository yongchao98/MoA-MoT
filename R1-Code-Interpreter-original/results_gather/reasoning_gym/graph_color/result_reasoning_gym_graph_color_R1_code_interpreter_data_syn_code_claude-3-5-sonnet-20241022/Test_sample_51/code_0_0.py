def is_valid_coloring(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(vertices, edges, available_colors):
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    
    # Initialize colors dictionary with None
    colors = {v: None for v in vertices}
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        current_vertex = vertices[vertex_index]
        for color in available_colors:
            if is_valid_coloring(graph, colors, current_vertex, color):
                colors[current_vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                colors[current_vertex] = None
        
        return False
    
    # Start backtracking from vertex 0
    backtrack(0)
    
    # Convert the result to the required format
    result = {str(k): v for k, v in colors.items()}
    print(f"<<<{result}>>>")

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (0, 8), (2, 4), (2, 5), (5, 8), (6, 9)]
available_colors = [1, 2, 3]

graph_coloring(vertices, edges, available_colors)