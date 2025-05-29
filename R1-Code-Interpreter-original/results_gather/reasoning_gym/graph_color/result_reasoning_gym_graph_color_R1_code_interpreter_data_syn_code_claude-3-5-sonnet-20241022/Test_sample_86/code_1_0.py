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

# Process vertices with connections first
connected_vertices = set()
for edge in edges:
    connected_vertices.add(edge[0])
    connected_vertices.add(edge[1])

# Color connected vertices first
for vertex in connected_vertices:
    if vertex not in coloring:
        used_colors = {coloring[adj] for adj in adj_list[vertex] if adj in coloring}
        for color in colors:
            if color not in used_colors:
                coloring[vertex] = color
                break

# Color remaining (isolated) vertices
for vertex in vertices:
    if vertex not in coloring:
        coloring[vertex] = 1  # isolated vertices can have any color, we use 1

# Convert to string format
result = {str(k): v for k, v in coloring.items()}
print(json.dumps(result))