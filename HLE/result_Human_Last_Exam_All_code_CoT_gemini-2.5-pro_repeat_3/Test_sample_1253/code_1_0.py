import collections

def solve():
    """
    Solves the path cover problem for the given examples.
    """
    # Example 1:
    # V = {S, 1, 2, 3, 4}, E = {(S, 1),  (S, 2), (1, 3), (2, 3), (2,4)} 
    # Expected output: 2
    graph1 = {
        'S': ['1', '2'],
        '1': ['S', '3'],
        '2': ['S', '3', '4'],
        '3': ['1', '2'],
        '4': ['2']
    }
    start_node1 = 'S'
    
    print("Analyzing Example 1:")
    find_min_path_cover(graph1, start_node1)
    
    print("\n" + "="*30 + "\n")

    # Example 2:
    # V = {S, 1, 2, 3, 4, 5, 6, 7, 8}, E = {(S, 1), (S, 2), (S, 3), (1, 4), (2, 4), (2, 5), (3, 5), (4, 6), (4, 7), (5, 8)
    # Expected output: 3
    graph2 = {
        'S': ['1', '2', '3'],
        '1': ['S', '4'],
        '2': ['S', '4', '5'],
        '3': ['S', '5'],
        '4': ['1', '2', '6', '7'],
        '5': ['2', '3', '8'],
        '6': ['4'],
        '7': ['4'],
        '8': ['5']
    }
    start_node2 = 'S'
    
    print("Analyzing Example 2:")
    find_min_path_cover(graph2, start_node2)

def find_min_path_cover(graph, start_node):
    """
    Calculates the minimum number of shortest paths from a start node to cover all vertices.
    This is equivalent to finding the size of the largest layer in a BFS traversal.
    
    Args:
        graph (dict): Adjacency list representation of the graph.
        start_node: The starting vertex for the paths.
    """
    if start_node not in graph:
        print("Start node not in graph.")
        return

    # Step 1 & 2: Perform BFS to find layers (distances)
    distances = {}
    queue = collections.deque([(start_node, 0)])
    distances[start_node] = 0

    while queue:
        current_node, dist = queue.popleft()
        
        for neighbor in graph.get(current_node, []):
            if neighbor not in distances:
                distances[neighbor] = dist + 1
                queue.append((neighbor, dist + 1))
    
    # Step 3: Count vertices in each layer
    layer_counts = collections.defaultdict(int)
    for node in distances:
        dist = distances[node]
        layer_counts[dist] += 1
        
    if not layer_counts:
        print("Graph is empty or start node has no connections.")
        return
        
    # Step 4: Find the maximum layer size
    max_layer_size = 0
    if layer_counts:
        max_layer_size = max(layer_counts.values())

    # Output the results as an equation
    print(f"Layer counts by distance from '{start_node}': {dict(sorted(layer_counts.items()))}")
    
    count_list = [str(c) for c in sorted(layer_counts.values(), reverse=True)]
    print(f"The size of the largest layer is: max({', '.join(count_list)}) = {max_layer_size}")
    
    print(f"\nThe size of the smallest set of paths is {max_layer_size}")


solve()