from collections import deque

def solve_path_cover(V, E, S):
    """
    Finds the size of the smallest set of S-rooted shortest paths
    that cover all vertices in the graph.

    Args:
        V (set): The set of vertices.
        E (list of tuples): The list of edges.
        S (any): The starting vertex.
    """
    if not V:
        print("0")
        return

    # Build adjacency list from edge list
    adj = {v: [] for v in V}
    for u, v in E:
        adj[u].append(v)
        adj[v].append(u)

    # Step 1: Separate vertices into layers by distance from S using BFS
    # dist will also serve as a 'visited' set
    dist = {v: -1 for v in V}
    layers = []

    q = deque([S])
    dist[S] = 0

    while True:
        # At the start of this loop, the queue contains all nodes
        # of the current layer.
        layer_nodes = list(q)
        if not layer_nodes:
            break
        layers.append(layer_nodes)
        
        # Process the current layer to find the next layer
        level_size = len(q)
        for _ in range(level_size):
            u = q.popleft()
            for v in adj[u]:
                if dist[v] == -1: # If v is not visited yet
                    dist[v] = dist[u] + 1
                    q.append(v)

    # Step 2: The result is the maximum size of any such layer.
    if not layers:
        max_layer_size = 0
    else:
        max_layer_size = max(len(layer) for layer in layers)

    # Print the layers for clarity
    print("Graph Layers by distance from S:")
    for i, layer in enumerate(layers):
        print(f"Layer {i} (distance {i}): {sorted(list(layer))}, size = {len(layer)}")
    
    print("\nThe size of the largest layer determines the minimum number of paths.")
    print(f"Result: {max_layer_size}")


if __name__ == '__main__':
    print("--- Example 1 ---")
    V1 = {'S', 1, 2, 3, 4}
    E1 = [('S', 1), ('S', 2), (1, 3), (2, 3), (2, 4)]
    S1 = 'S'
    solve_path_cover(V1, E1, S1)
    
    print("\n--- Example 2 ---")
    V2 = {'S', 1, 2, 3, 4, 5, 6, 7, 8}
    E2 = [('S', 1), ('S', 2), ('S', 3), (1, 4), (2, 4), (2, 5), (3, 5), 
          (4, 6), (4, 7), (5, 8)]
    S2 = 'S'
    solve_path_cover(V2, E2, S2)
