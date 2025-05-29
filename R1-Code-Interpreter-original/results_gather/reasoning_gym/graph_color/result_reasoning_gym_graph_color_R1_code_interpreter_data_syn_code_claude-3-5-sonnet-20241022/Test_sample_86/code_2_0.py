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

# Find connected components
def find_component(vertex, visited):
    component = []
    stack = [vertex]
    while stack:
        v = stack.pop()
        if v not in visited:
            visited.add(v)
            component.append(v)
            stack.extend(adj for adj in adj_list[v] if adj not in visited)
    return component

# Color a component
def color_component(component):
    for vertex in component:
        used_colors = {coloring[adj] for adj in adj_list[vertex] if adj in coloring}
        for color in colors:
            if color not in used_colors:
                coloring[vertex] = color
                break

# Process vertices by components
visited = set()
for vertex in vertices:
    if vertex not in visited:
        component = find_component(vertex, visited)
        if len(component) > 1:  # If component has connected vertices
            color_component(component)
        else:  # Isolated vertex
            coloring[vertex] = 1

# Convert to string format
result = {str(k): v for k, v in coloring.items()}
print(json.dumps(result))