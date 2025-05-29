import json

def is_valid_coloring(graph, colors):
    # Check if the current coloring is valid
    for v1, neighbors in graph.items():
        for v2 in neighbors:
            if v1 in colors and v2 in colors and colors[v1] == colors[v2]:
                return False
    return True

# Create graph representation (adjacency list)
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (1, 6), (4, 7), (7, 8)]

# Create adjacency list
graph = {v: [] for v in vertices}
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Initialize with vertices that have edges
connected_vertices = set()
for edge in edges:
    connected_vertices.add(edge[0])
    connected_vertices.add(edge[1])

# Color connected vertices first
colors = {}
for vertex in connected_vertices:
    used_colors = {colors[neighbor] for neighbor in graph[vertex] if neighbor in colors}
    available_colors = [c for c in [1, 2, 3] if c not in used_colors]
    colors[vertex] = available_colors[0]

# Color remaining vertices
remaining_vertices = set(vertices) - connected_vertices
for vertex in remaining_vertices:
    colors[vertex] = 1  # Can use any color for disconnected vertices

# Verify the coloring
assert is_valid_coloring(graph, colors)

# Convert to string format required
result = {str(k): v for k, v in sorted(colors.items())}
print(f"<<<{json.dumps(result)}>>>")