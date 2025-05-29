import json

# Create graph representation
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (0, 6), (5, 7)]
colors = [1, 2, 3]

# Create adjacency list
adj_list = {v: [] for v in vertices}
for v1, v2 in edges:
    adj_list[v1].append(v2)
    adj_list[v2].append(v1)

# Initialize coloring
coloring = {}

# Color vertices
for vertex in vertices:
    # For isolated vertices, just use color 1
    if not adj_list[vertex]:
        coloring[vertex] = 1
        continue
    
    # For connected vertices, use first available color
    used_colors = {coloring[adj] for adj in adj_list[vertex] if adj in coloring}
    for color in colors:
        if color not in used_colors:
            coloring[vertex] = color
            break

# Convert to string format
result = {str(k): v for k, v in coloring.items()}
print(json.dumps(result))