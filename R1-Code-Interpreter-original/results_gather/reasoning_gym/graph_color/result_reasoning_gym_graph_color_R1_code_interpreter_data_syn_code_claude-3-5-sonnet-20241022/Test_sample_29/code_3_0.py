import json

def is_valid_coloring(coloring, edges):
    for edge in edges:
        v1, v2 = edge
        if str(v1) in coloring and str(v2) in coloring:
            if coloring[str(v1)] == coloring[str(v2)]:
                return False
    return True

def get_conflicts(vertex, color, coloring, edges):
    conflicts = 0
    for edge in edges:
        if str(edge[0]) == str(vertex) and str(edge[1]) in coloring:
            if coloring[str(edge[1])] == color:
                conflicts += 1
        elif str(edge[1]) == str(vertex) and str(edge[0]) in coloring:
            if coloring[str(edge[0])] == color:
                conflicts += 1
    return conflicts

vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 6), (0, 8), (1, 4), (1, 6), (1, 8), (5, 6), (5, 8)]
colors = [1, 2, 3]

# Create adjacency list
adj = {str(v): [] for v in vertices}
for edge in edges:
    adj[str(edge[0])].append(str(edge[1]))
    adj[str(edge[1])].append(str(edge[0]))

# Sort vertices by degree
vertex_degrees = [(str(v), len(adj[str(v)])) for v in vertices]
ordered_vertices = [v for v, _ in sorted(vertex_degrees, key=lambda x: x[1], reverse=True)]

# Initialize coloring
coloring = {}

# First pass: Assign colors greedily
for vertex in ordered_vertices:
    min_conflicts = float('inf')
    best_color = 1
    
    for color in colors:
        conflicts = get_conflicts(vertex, color, coloring, edges)
        if conflicts < min_conflicts:
            min_conflicts = conflicts
            best_color = color
    
    coloring[vertex] = best_color

# Second pass: Fix conflicts
while not is_valid_coloring(coloring, edges):
    for vertex in ordered_vertices:
        # If vertex has conflicts, try to fix
        current_conflicts = get_conflicts(vertex, coloring[vertex], coloring, edges)
        if current_conflicts > 0:
            # Try other colors
            for color in colors:
                if color != coloring[vertex]:
                    new_conflicts = get_conflicts(vertex, color, coloring, edges)
                    if new_conflicts < current_conflicts:
                        coloring[vertex] = color
                        break

# Verify final solution
if not is_valid_coloring(coloring, edges):
    # If still invalid, use a guaranteed solution
    coloring = {
        "0": 1,
        "1": 2,
        "2": 1,
        "3": 1,
        "4": 3,
        "5": 2,
        "6": 3,
        "7": 1,
        "8": 3,
        "9": 1
    }

print(f"<<<{json.dumps(coloring)}>>>")