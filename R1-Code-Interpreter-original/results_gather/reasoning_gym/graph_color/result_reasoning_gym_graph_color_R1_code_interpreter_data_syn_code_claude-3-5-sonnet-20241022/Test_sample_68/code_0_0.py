import json

# Create graph representation
vertices = list(range(10))
edges = [(5, 6), (7, 9)]
colors = [1, 2, 3]

# Initialize coloring
coloring = {}

# Create adjacency list for easier processing
adj_list = {v: [] for v in vertices}
for v1, v2 in edges:
    adj_list[v1].append(v2)
    adj_list[v2].append(v1)

# Color vertices
for v in vertices:
    # For vertices with no neighbors, just use color 1
    if not adj_list[v]:
        coloring[v] = 1
        continue
    
    # For vertices with neighbors, find first available color
    neighbor_colors = {coloring[n] for n in adj_list[v] if n in coloring}
    for color in colors:
        if color not in neighbor_colors:
            coloring[v] = color
            break

# Convert to string format
result = {str(k): v for k, v in coloring.items()}
print("<<<" + json.dumps(result) + ">>>")