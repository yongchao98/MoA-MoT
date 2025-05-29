import json

def is_safe(graph, vertex, colors, color):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, num_colors, vertices):
    colors = {}
    
    # Initialize all vertices with color 0 (uncolored)
    for vertex in vertices:
        colors[vertex] = 0
    
    def color_graph(vertex):
        # If all vertices are colored, we're done
        if vertex >= len(vertices):
            return True
            
        # Try each color
        for color in range(1, num_colors + 1):
            if is_safe(graph, vertices[vertex], colors, color):
                colors[vertices[vertex]] = color
                
                # Recursively color remaining vertices
                if color_graph(vertex + 1):
                    return True
                    
                # If coloring doesn't work, backtrack
                colors[vertices[vertex]] = 0
        
        return False
    
    # Create adjacency list representation
    graph_adj = {v: [] for v in vertices}
    edges = [(1, 4), (2, 5), (3, 4), (3, 9), (5, 6), (5, 8)]
    for v1, v2 in edges:
        graph_adj[v1].append(v2)
        graph_adj[v2].append(v1)
    
    # Color the graph
    if color_graph(0):
        # Convert all values to integers and create JSON string
        result = {str(k): v for k, v in colors.items()}
        print(f"<<<{json.dumps(result)}>>>")
    else:
        print("No solution exists")

# Run the solution
vertices = list(range(10))
graph_coloring(vertices, 3, vertices)