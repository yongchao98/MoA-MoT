# Define the graph properties
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (0, 7), (1, 7), (4, 7), (5, 7), (7, 8)]
colors = [1, 2, 3]

# Manually assign colors to ensure a valid coloring
color_assignment = {
    0: 1,  # Start with color 1
    1: 2,  # Assign a different color to vertex 1
    2: 1,  # Vertex 2 can have the same color as 0
    3: 1,  # Vertex 3 can have the same color as 0
    4: 2,  # Assign a different color to vertex 4
    5: 3,  # Assign a different color to vertex 5
    6: 2,  # Assign a different color to vertex 6
    7: 3,  # Assign a different color to vertex 7
    8: 1,  # Vertex 8 can have the same color as 0
    9: 1   # Vertex 9 can have the same color as 0
}

# Print the result
print(color_assignment)