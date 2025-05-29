def is_valid_coloring(graph, colors):
    # Check if adjacent vertices have different colors
    for v1, v2 in graph:
        if str(v1) in colors and str(v2) in colors:
            if colors[str(v1)] == colors[str(v2)]:
                return False
    return True

# Create graph
edges = [(1, 7), (2, 4), (6, 7)]
vertices = list(range(10))
possible_colors = [1, 2, 3]

# Initialize coloring
coloring = {}

# For isolated vertices, use color 1
for v in vertices:
    coloring[str(v)] = 1

# Handle connected components
for v1, v2 in edges:
    # If vertices are already colored differently, continue
    if coloring[str(v1)] != coloring[str(v2)]:
        continue
    # If they have same color, change second vertex color
    coloring[str(v2)] = 2

# Verify the coloring is valid
if is_valid_coloring(edges, coloring):
    print(str(coloring).replace("'", '"'))