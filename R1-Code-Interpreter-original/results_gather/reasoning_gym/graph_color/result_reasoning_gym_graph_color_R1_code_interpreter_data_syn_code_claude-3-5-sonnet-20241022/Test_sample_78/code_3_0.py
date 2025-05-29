import json

def is_safe(graph, vertex, color, colors):
    # Check if it's safe to color vertex with given color
    for neighbor in graph[vertex]:
        if neighbor in colors and colors[neighbor] == color:
            return False
    return True

def graph_coloring(graph, vertices, m_colors, vertex_index, colors):
    # Base case: if all vertices are colored
    if vertex_index == len(vertices):
        return True
    
    vertex = vertices[vertex_index]
    # Try all colors
    for color in range(1, m_colors + 1):
        # Check if assignment of color is safe
        if is_safe(graph, vertex, color, colors):
            colors[vertex] = color
            # Recur to assign colors to remaining vertices
            if graph_coloring(graph, vertices, m_colors, vertex_index + 1, colors):
                return True
            # If color doesn't lead to solution, remove it
            colors[vertex] = 0
    
    return False

# Create graph representation
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(1, 4), (1, 6), (4, 7), (7, 8)]
m_colors = 3

# Create adjacency list
graph = {v: [] for v in vertices}
for v1, v2 in edges:
    graph[v1].append(v2)
    graph[v2].append(v1)

# Initialize colors dictionary
colors = {v: 0 for v in vertices}

# Find coloring using backtracking
success = graph_coloring(graph, vertices, m_colors, 0, colors)

# Verify the coloring
valid = True
for v1, v2 in edges:
    if colors[v1] == colors[v2]:
        valid = False
        break

assert valid and success, "Invalid coloring found"

# Convert to string format required
result = {str(k): v for k, v in sorted(colors.items())}
print(f"<<<{json.dumps(result)}>>>")