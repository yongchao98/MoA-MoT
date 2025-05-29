# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 7), (2, 4), (2, 8), (2, 9), (4, 6), (5, 8), (5, 9), (7, 9)]
colors = [1, 2, 3]

# Initialize the color assignment dictionary
color_assignment = {vertex: None for vertex in vertices}

# Create an adjacency list for the graph
adjacency_list = {vertex: [] for vertex in vertices}
for edge in edges:
    adjacency_list[edge[0]].append(edge[1])
    adjacency_list[edge[1]].append(edge[0])

# Function to check if a color can be assigned to a vertex
def is_safe(vertex, color):
    for neighbor in adjacency_list[vertex]:
        if color_assignment[neighbor] == color:
            return False
    return True

# Function to solve the graph coloring problem using backtracking
def graph_coloring(vertex):
    if vertex == len(vertices):
        return True
    
    for color in colors:
        if is_safe(vertex, color):
            color_assignment[vertex] = color
            if graph_coloring(vertex + 1):
                return True
            color_assignment[vertex] = None
    
    return False

# Start the coloring process
graph_coloring(0)

# Print the color assignment
print(color_assignment)