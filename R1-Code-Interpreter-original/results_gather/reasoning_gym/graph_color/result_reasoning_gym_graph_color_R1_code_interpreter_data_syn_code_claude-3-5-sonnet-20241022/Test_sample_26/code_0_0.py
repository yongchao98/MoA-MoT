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

# Color vertices
for vertex in vertices:
    # Get colors of adjacent vertices
    adj_colors = {coloring[adj] for adj in adj_list[vertex]}
    # Choose the first available color
    for color in colors:
        if color not in adj_colors:
            coloring[vertex] = color
            break

# Convert to string format and print
coloring_str = {str(k): v for k, v in coloring.items()}
print(json.dumps(coloring_str))