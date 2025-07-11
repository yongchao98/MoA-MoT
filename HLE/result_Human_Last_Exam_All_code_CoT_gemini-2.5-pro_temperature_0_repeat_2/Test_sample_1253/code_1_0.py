import collections

def solve():
    """
    Solves the minimum path cover problem for a given graph.
    The solution is based on finding the size of the largest layer in the 
    shortest-path graph from the source S, as explained by Dilworth's Theorem
    for graded posets.
    """
    # Example 2 Graph Definition
    vertices = {'S', '1', '2', '3', '4', '5', '6', '7', '8'}
    edges = [('S', '1'), ('S', '2'), ('S', '3'), ('1', '4'), ('2', '4'), 
             ('2', '5'), ('3', '5'), ('4', '6'), ('4', '7'), ('5', '8')]
    start_node = 'S'

    # Build adjacency list for the graph
    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # Step 1 & 2: Run BFS to find distances and group nodes into layers
    q = collections.deque([(start_node, 0)])
    distances = {start_node: 0}
    layers = collections.defaultdict(list)
    layers[0].append(start_node)
    
    max_dist = 0

    while q:
        u, dist = q.popleft()
        max_dist = max(max_dist, dist)

        for v in adj[u]:
            if v not in distances:
                distances[v] = dist + 1
                layers[dist + 1].append(v)
                q.append((v, dist + 1))

    # Step 3: Find the maximum size of any layer
    layer_sizes = []
    for i in range(max_dist + 1):
        # Sort nodes within a layer for consistent output
        layers[i].sort()
        layer_sizes.append(len(layers[i]))
        print(f"Layer {i} (distance from S): {layers[i]}, size = {len(layers[i])}")

    max_size = max(layer_sizes)
    
    # Print the final equation and result
    layer_sizes_str = ", ".join(map(str, layer_sizes))
    print(f"\nThe minimum number of paths required is the size of the largest layer.")
    print(f"Result = max({layer_sizes_str}) = {max_size}")

solve()