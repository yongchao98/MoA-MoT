import json

# Create graph representation
vertices = list(range(10))
edges = [(0, 4), (1, 8), (2, 8)]
colors = [1, 2, 3]

# Initialize all vertices with color 1
coloring = {v: 1 for v in vertices}

# Create adjacency list
adj_list = {v: [] for v in vertices}
for u, v in edges:
    adj_list[u].append(v)
    adj_list[v].append(u)

# First, color vertex 8 since it has the most constraints (connected to both 1 and 2)
coloring[8] = 1

# Color vertices 1 and 2 (must be different from vertex 8's color)
coloring[1] = 2
coloring[2] = 2  # Can be same as vertex 1 since they're not connected

# Color vertices 0 and 4
coloring[0] = 1
coloring[4] = 2  # Must be different from vertex 0

# Convert to string format and print
coloring_str = {str(k): v for k, v in coloring.items()}
print(json.dumps(coloring_str))