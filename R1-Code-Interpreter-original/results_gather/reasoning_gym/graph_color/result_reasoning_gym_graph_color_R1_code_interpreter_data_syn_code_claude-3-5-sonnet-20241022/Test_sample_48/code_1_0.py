import json

def create_graph(vertices, edges):
    graph = {v: [] for v in vertices}
    for edge in edges:
        graph[edge[0]].append(edge[1])
        graph[edge[1]].append(edge[0])
    return graph

def is_valid_coloring(graph, vertex, color, colors):
    # Check if any adjacent vertex has the same color
    for neighbor in graph[vertex]:
        if colors.get(neighbor) == color:
            return False
    return True

def find_coloring(vertices, edges, available_colors):
    graph = create_graph(vertices, edges)
    colors = {}
    
    def backtrack(vertex_idx):
        if vertex_idx == len(vertices):
            return True
            
        current_vertex = vertices[vertex_idx]
        
        # Try each color
        for color in available_colors:
            if is_valid_coloring(graph, current_vertex, color, colors):
                colors[current_vertex] = color
                if backtrack(vertex_idx + 1):
                    return True
                colors[current_vertex] = None  # Backtrack
        
        return False
    
    # Start the coloring process
    if backtrack(0):
        result = {str(k): v for k, v in colors.items()}
        print(json.dumps(result))
    else:
        print("No valid coloring found")

# Problem parameters
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 7), (1, 8), (4, 7), (5, 9), (8, 9)]
available_colors = [1, 2, 3]

# Find solution
find_coloring(vertices, edges, available_colors)