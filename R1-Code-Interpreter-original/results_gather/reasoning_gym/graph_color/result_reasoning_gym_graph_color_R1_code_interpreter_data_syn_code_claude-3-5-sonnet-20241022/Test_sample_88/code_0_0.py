def is_valid_coloring(graph, colors):
    # Check if adjacent vertices have different colors
    for v1, v2 in graph:
        if v1 in colors and v2 in colors and colors[v1] == colors[v2]:
            return False
    return True

def find_coloring():
    vertices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    edges = [(0, 6), (2, 9), (4, 9), (5, 6)]
    possible_colors = [1, 2, 3]
    
    # Initialize all vertices with color 1
    coloring = {v: 1 for v in vertices}
    
    # Create adjacency list
    adj_list = {}
    for v in vertices:
        adj_list[v] = []
    for v1, v2 in edges:
        adj_list[v1].append(v2)
        adj_list[v2].append(v1)
    
    # For each vertex with edges, ensure valid coloring
    for v in vertices:
        if adj_list[v]:
            # Get colors of neighbors
            neighbor_colors = {coloring[n] for n in adj_list[v]}
            # Find first available color
            for color in possible_colors:
                if color not in neighbor_colors:
                    coloring[v] = color
                    break
    
    # Convert to JSON string format
    result = str(coloring).replace("'", '"')
    print(f"<<<{result}>>>")

find_coloring()