import json

def get_valid_color(graph, vertex, colors, available_colors):
    # Get colors used by neighbors
    neighbor_colors = set()
    for neighbor in graph[vertex]:
        if colors[neighbor] != 0:  # if neighbor is colored
            neighbor_colors.add(colors[neighbor])
    
    # Return first available color that's not used by neighbors
    for color in available_colors:
        if color not in neighbor_colors:
            return color
    return None

# Create graph representation
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 4), (1, 2), (1, 3), (1, 9), (3, 5), (4, 8), (7, 9)]
available_colors = [1, 2, 3]

# Create adjacency list
graph = {v: [] for v in vertices}
for edge in edges:
    graph[edge[0]].append(edge[1])
    graph[edge[1]].append(edge[0])

# Initialize colors
colors = {v: 0 for v in vertices}

# Process vertices in order of degree (number of neighbors)
vertex_degrees = [(v, len(graph[v])) for v in vertices]
vertex_degrees.sort(key=lambda x: x[1], reverse=True)

# Color vertices
for vertex, _ in vertex_degrees:
    color = get_valid_color(graph, vertex, colors, available_colors)
    colors[vertex] = color if color is not None else 1

# Verify solution
valid = True
for v in vertices:
    for neighbor in graph[v]:
        if colors[v] == colors[neighbor]:
            valid = False
            break
    if not valid:
        break

# If solution is not valid, try a different approach
if not valid:
    # Reset colors
    colors = {v: 0 for v in vertices}
    
    # Start with vertex 1 (which has most connections)
    colors[1] = 1
    for neighbor in graph[1]:
        colors[neighbor] = 2  # Give all neighbors of vertex 1 color 2
    
    # Color remaining vertices
    for v in vertices:
        if colors[v] == 0:
            colors[v] = get_valid_color(graph, v, colors, available_colors) or 3

# Convert to string keys for JSON
solution = {str(k): v for k, v in colors.items()}

print(f"<<<{json.dumps(solution)}>>>")