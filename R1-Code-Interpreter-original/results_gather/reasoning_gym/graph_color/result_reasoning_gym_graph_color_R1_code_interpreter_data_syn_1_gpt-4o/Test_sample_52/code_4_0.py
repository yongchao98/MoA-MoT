def is_valid(vertex, color, color_assignment, adjacency_list):
    for neighbor in adjacency_list[vertex]:
        if color_assignment[neighbor] == color:
            return False
    return True

def graph_coloring_backtrack(vertices, edges, colors):
    # Initialize the color assignment for each vertex
    color_assignment = {v: None for v in vertices}
    
    # Create an adjacency list for the graph
    adjacency_list = {v: [] for v in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    def backtrack(vertex_index):
        if vertex_index == len(vertices):
            return True
        
        vertex = vertices[vertex_index]
        for color in colors:
            if is_valid(vertex, color, color_assignment, adjacency_list):
                color_assignment[vertex] = color
                if backtrack(vertex_index + 1):
                    return True
                color_assignment[vertex] = None
        
        return False
    
    backtrack(0)
    return color_assignment

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 2), (0, 3), (2, 5), (3, 5), (4, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring_backtrack(vertices, edges, colors)
print(color_assignment)