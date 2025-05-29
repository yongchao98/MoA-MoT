def is_valid_color(graph, colors, vertex, color):
    # Check if any adjacent vertex has the same color
    for adj_vertex in range(len(graph)):
        if graph[vertex][adj_vertex] == 1 and colors[adj_vertex] == color:
            return False
    return True

def graph_coloring(graph, num_colors, vertices, colors, vertex):
    # Base case: if all vertices are colored
    if vertex == len(vertices):
        return True
    
    # Try each color for the current vertex
    for color in range(1, num_colors + 1):
        if is_valid_color(graph, colors, vertex, color):
            colors[vertex] = color
            if graph_coloring(graph, num_colors, vertices, colors, vertex + 1):
                return True
            colors[vertex] = 0
    
    return False

# Create adjacency matrix from edges
vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
edges = [(0, 5), (0, 6), (0, 8), (1, 3), (3, 9), (4, 5), (4, 6), (6, 8), (7, 8)]
num_vertices = len(vertices)
graph = [[0] * num_vertices for _ in range(num_vertices)]

for edge in edges:
    graph[edge[0]][edge[1]] = 1
    graph[edge[1]][edge[0]] = 1

# Initialize colors array
colors = [0] * num_vertices

# Find coloring
graph_coloring(graph, 3, vertices, colors, 0)

# Convert result to required format
result = {str(i): colors[i] for i in range(num_vertices)}
print(f"<<<{result}>>>")