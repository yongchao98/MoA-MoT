def is_safe(vertex, color, color_assignment, adjacency_list):
    # Check if the color can be assigned to the vertex
    for neighbor in adjacency_list[vertex]:
        if color_assignment[neighbor] == color:
            return False
    return True

def graph_coloring_util(vertices, colors, color_assignment, adjacency_list, vertex_index):
    # If all vertices are assigned a color, return True
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    
    # Try different colors for the vertex
    for color in colors:
        if is_safe(vertex, color, color_assignment, adjacency_list):
            color_assignment[vertex] = color
            
            # Recur to assign colors to the rest of the vertices
            if graph_coloring_util(vertices, colors, color_assignment, adjacency_list, vertex_index + 1):
                return True
            
            # If assigning color doesn't lead to a solution, remove it
            color_assignment[vertex] = None
    
    return False

def graph_coloring(vertices, edges, colors):
    # Initialize the color assignment dictionary
    color_assignment = {vertex: None for vertex in vertices}
    
    # Create an adjacency list from the edges
    adjacency_list = {vertex: [] for vertex in vertices}
    for edge in edges:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    
    # Start the coloring process
    if graph_coloring_util(vertices, colors, color_assignment, adjacency_list, 0):
        return color_assignment
    else:
        return None

# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 2), (1, 8), (2, 3), (3, 4), (8, 9)]
colors = [1, 2, 3]

# Get the color assignment
color_assignment = graph_coloring(vertices, edges, colors)
print(color_assignment)