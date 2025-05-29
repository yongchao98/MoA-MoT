def is_safe(graph, vertex, color, coloring):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if neighbor in coloring and coloring[neighbor] == color:
            return False
    return True

def graph_coloring(graph, colors, vertices):
    # Initialize empty coloring dictionary
    coloring = {}
    
    # Create adjacency list representation
    adj_list = {v: [] for v in vertices}
    for edge in [(0, 1), (0, 6), (5, 9), (6, 8), (8, 9)]:
        adj_list[edge[0]].append(edge[1])
        adj_list[edge[1]].append(edge[0])
    
    def backtrack(vertex_idx):
        if vertex_idx == len(vertices):
            return True
        
        current_vertex = vertices[vertex_idx]
        
        # Try each color
        for color in colors:
            if is_safe(adj_list, current_vertex, color, coloring):
                coloring[current_vertex] = color
                if backtrack(vertex_idx + 1):
                    return True
                coloring.pop(current_vertex, None)
        
        return False
    
    # Start backtracking from first vertex
    backtrack(0)
    
    # Convert all values to strings for JSON compatibility
    result = {str(k): v for k, v in coloring.items()}
    
    print(f"<<<{result}>>>")

# Input data
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
colors = [1, 2, 3]

graph_coloring({}, colors, vertices)