def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in colors:
            if is_safe(graph, vertex, color, coloring):
                coloring[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
        return False
    
    # Initialize graph as adjacency list
    graph_adj = {v: [] for v in vertices}
    for edge in [(0, 9), (2, 3), (2, 6), (3, 8), (5, 7)]:
        graph_adj[edge[0]].append(edge[1])
        graph_adj[edge[1]].append(edge[0])
    
    # Initialize coloring dictionary
    coloring = {}
    
    # Start backtracking from vertex 0
    backtrack(0)
    
    # Convert result to JSON string
    result = {str(k): v for k, v in coloring.items()}
    print("<<<" + str(result) + ">>>")

# Run the solution
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
colors = [1, 2, 3]
graph_coloring({v: [] for v in vertices}, colors, vertices)