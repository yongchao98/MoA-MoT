import collections

def solve():
    """
    Solves the path cover problem for two example graphs.
    """

    # --- Example 1 ---
    # V = {S, 1, 2, 3, 4}, E = {(S, 1),  (S, 2), (1, 3), (2, 3), (2,4)}
    # Expected output: 2
    edges1 = [
        ('S', '1'), ('S', '2'),
        ('1', '3'),
        ('2', '3'), ('2', '4')
    ]
    start_node1 = 'S'
    print("--- Example 1 ---")
    find_min_path_cover_size(edges1, start_node1)
    print("\n" + "="*20 + "\n")

    # --- Example 2 ---
    # V = {S, 1..8}, E = {(S,1),(S,2),(S,3), (1,4),(2,4),(2,5),(3,5),(4,6),(4,7),(5,8)}
    # Expected output: 3
    edges2 = [
        ('S', '1'), ('S', '2'), ('S', '3'),
        ('1', '4'),
        ('2', '4'), ('2', '5'),
        ('3', '5'),
        ('4', '6'), ('4', '7'),
        ('5', '8')
    ]
    start_node2 = 'S'
    print("--- Example 2 ---")
    find_min_path_cover_size(edges2, start_node2)


def find_min_path_cover_size(edges, start_node):
    """
    Calculates the size of the smallest set of paths based on layer sizes.

    Args:
        edges: A list of tuples representing the graph edges.
        start_node: The starting vertex S.
    """
    if not edges:
        print("Result: 0")
        return

    # Build adjacency list from edges
    adj = collections.defaultdict(list)
    nodes = set()
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        nodes.add(u)
        nodes.add(v)

    if start_node not in nodes:
        print(f"Start node {start_node} not in graph.")
        return

    # Step 1 & 2: Run BFS and find layers
    queue = collections.deque([(start_node, 0)])
    visited = {start_node}
    # Store nodes by their layer (distance from start_node)
    layers = collections.defaultdict(list)
    layers[0].append(start_node)
    
    max_layer_size = 0

    while queue:
        current_node, distance = queue.popleft()

        for neighbor in adj[current_node]:
            if neighbor not in visited:
                visited.add(neighbor)
                new_distance = distance + 1
                layers[new_distance].append(neighbor)
                queue.append((neighbor, new_distance))

    # Step 3: Find the maximum size of any layer (excluding layer 0)
    print("Layers by distance from S:")
    # Start from layer 1 since layer 0 is just the start node
    if len(layers) > 1:
        # Find the max size of layers 1 and onwards
        max_layer_size = max(len(layer_nodes) for i, layer_nodes in layers.items() if i > 0)
        
    for i in sorted(layers.keys()):
      print(f"Layer {i} (size {len(layers[i])}): {sorted(layers[i])}")


    print(f"\nThe maximum size of any layer (excluding layer 0) is {max_layer_size}.")
    print(f"Final Answer: The size of the smallest set of paths is {max_layer_size}")

solve()