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
    
    # Initialize colors dictionary
    colors = {v: None for v in vertices}
    
    def backtrack(vertex):
        # If we've colored all vertices, we're done
        if vertex == len(vertices):
            return True
        
        # Try each color for current vertex
        for color in available_colors:
            if is_valid_coloring(graph, colors, vertices[vertex], color):
                colors[vertices[vertex]] = color
                if backtrack(vertex + 1):
                    return True
                colors[vertices[vertex]] = None
        
        return False
    
    # Start backtracking from first vertex
    backtrack(0)
    
    # Convert result to string format
    result = {str(k): v for k, v in colors.items()}
    print("<<<" + str(result) + ">>>")

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (1, 4), (2, 6), (2, 7), (3, 7)]
colors = [1, 2, 3]

graph_coloring(vertices, edges, colors)