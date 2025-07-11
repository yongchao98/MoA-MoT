import collections

def solve_path_cover():
    """
    Calculates the size of the smallest set of paths to cover the graph
    based on algorithm A: finding the maximum layer size from the source S.
    """
    # Example 2 graph definition
    # Using strings for node names to handle 'S' gracefully
    edges = [
        ('S', '1'), ('S', '2'), ('S', '3'), 
        ('1', '4'), 
        ('2', '4'), ('2', '5'), 
        ('3', '5'), 
        ('4', '6'), ('4', '7'), 
        ('5', '8')
    ]
    source = 'S'

    # Build adjacency list for the undirected graph
    adj = collections.defaultdict(list)
    nodes = set()
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        nodes.add(u)
        nodes.add(v)

    # Step 1: Separate vertices into layers by distance from S using BFS
    # A queue for BFS, storing (node, distance)
    queue = collections.deque([(source, 0)])
    # A set to keep track of visited nodes
    visited = {source}
    # A dictionary to store layers: {distance: [nodes]}
    layers = collections.defaultdict(list)
    layers[0].append(source)

    while queue:
        current_node, distance = queue.popleft()
        
        for neighbor in adj[current_node]:
            if neighbor not in visited:
                visited.add(neighbor)
                new_distance = distance + 1
                layers[new_distance].append(neighbor)
                queue.append((neighbor, new_distance))
    
    # Check if all nodes were reached
    if len(visited) != len(nodes):
        print("Warning: Not all nodes are reachable from the source.")

    # Step 2: The result is the maximum size of any layer
    max_size = 0
    if not layers:
        print("Result: 0")
        return

    print("Layer sizes calculated from S:")
    layer_sizes = []
    # Iterate through sorted layer distances for clear output
    for dist in sorted(layers.keys()):
        size = len(layers[dist])
        layer_sizes.append(str(size))
        max_size = max(max_size, size)
    
    # As requested, output the numbers in the "final equation"
    # This shows the sizes of the layers found.
    print(f"max({', '.join(layer_sizes)}) = {max_size}")

    # The final answer
    print(f"\nThe size of the smallest set of paths is: {max_size}")


if __name__ == "__main__":
    solve_path_cover()