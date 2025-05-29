import json

# Create graph representation
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (1, 6), (4, 7), (7, 8)]

# Create adjacency list
graph = {v: [] for v in vertices}
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Initialize colors dictionary
colors = {v: 1 for v in vertices}  # Default all vertices to color 1

# Explicitly color the connected components
# Path: 1 - 4 - 7 - 8
colors[1] = 1
colors[4] = 2
colors[7] = 3
colors[8] = 2

# Handle vertex 6 which is connected to 1
colors[6] = 3

# Verify the coloring
valid = True
for v1, v2 in edges:
    if colors[v1] == colors[v2]:
        valid = False
        break

assert valid, "Invalid coloring found"

# Convert to string format required
result = {str(k): v for k, v in sorted(colors.items())}
print(f"<<<{json.dumps(result)}>>>")