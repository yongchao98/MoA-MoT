import json

# Create graph representation
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (0, 6), (5, 7)]
colors = [1, 2, 3]

# Initialize all vertices with color 1 (default color for isolated vertices)
coloring = {v: 1 for v in vertices}

# Handle first component (0, 4, 6) explicitly
coloring[0] = 1  # Assign first color to central vertex
coloring[4] = 2  # Assign second color to first neighbor
coloring[6] = 3  # Assign third color to second neighbor

# Handle second component (5, 7) explicitly
coloring[5] = 1  # Assign first color
coloring[7] = 2  # Assign second color

# Convert to string format
result = {str(k): v for k, v in coloring.items()}
print(json.dumps(result))