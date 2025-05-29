# Define the graph
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 9), (2, 4), (3, 9), (5, 6), (6, 7)]
colors = [1, 2, 3]

# Initialize the color assignment
color_assignment = {vertex: None for vertex in vertices}

# Create an adjacency list
adjacency_list = {vertex: [] for vertex in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Function to check if the current color assignment is valid
def is_valid(vertex, color):
    for neighbor in adjacency_list[vertex]:
        if color_assignment[neighbor] == color:
            return False
    return True

# Backtracking function to assign colors
def assign_colors(vertex_index):
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    for color in colors:
        if is_valid(vertex, color):
            color_assignment[vertex] = color
            if assign_colors(vertex_index + 1):
                return True
            color_assignment[vertex] = None
    
    return False

# Start the coloring process
assign_colors(0)

print(color_assignment)