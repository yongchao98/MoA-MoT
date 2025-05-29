import json

def is_safe(graph, vertex, colors, color):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring():
    # Initialize the graph
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(1, 4), (2, 5), (3, 4), (3, 9), (5, 6), (5, 8)]
    colors = {}
    num_colors = 3
    
    # Create adjacency list representation
    graph = {v: [] for v in vertices}
    for v1, v2 in edges:
        graph[v1].append(v2)
        graph[v2].append(v1)
    
    def color_vertex(vertex_idx):
        if vertex_idx == len(vertices):
            return True
            
        current_vertex = vertices[vertex_idx]
        for color in range(1, num_colors + 1):
            if is_safe(graph, current_vertex, colors, color):
                colors[current_vertex] = color
                if color_vertex(vertex_idx + 1):
                    return True
                colors.pop(current_vertex, None)
        
        return False
    
    # Try to color the graph
    if color_vertex(0):
        # Convert the solution to the required format
        result = {str(v): colors.get(v, 1) for v in vertices}  # Default color 1 for disconnected vertices
        print(f"<<<{json.dumps(result)}>>>")
    else:
        print("No solution exists")

# Run the solution
graph_coloring()