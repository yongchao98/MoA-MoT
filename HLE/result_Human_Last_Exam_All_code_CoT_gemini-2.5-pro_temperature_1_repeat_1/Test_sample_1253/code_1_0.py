import collections

def solve_path_cover(V_count, edges, S_node):
    """
    Calculates the size of the smallest set of shortest paths from S to cover all vertices.

    The algorithm works as follows:
    1.  Perform a Breadth-First Search (BFS) starting from the source node S to
        determine the shortest distance of every other node from S. This groups
        the nodes into "layers", where all nodes in a layer are equidistant from S.
    2.  Count the number of nodes in each layer.
    3.  The size of the minimum path cover is the size of the largest layer. This is
        because any single shortest path can only visit at most one node from any
        given layer. Therefore, to cover the largest layer, we need a number of
        paths equal to its size. This number is also proven to be sufficient.
    """
    if V_count == 0:
        print("Graph is empty.")
        return

    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # Step 1: Perform BFS to find layers
    q = collections.deque([(S_node, 0)])
    distances = {S_node: 0}
    layer_counts = collections.defaultdict(int)
    layer_counts[0] = 1

    while q:
        u, dist = q.popleft()
        for v in adj[u]:
            if v not in distances:
                distances[v] = dist + 1
                layer_counts[dist + 1] += 1
                q.append((v, dist + 1))
    
    if not layer_counts:
        # Handle the case where S is an isolated vertex
        if V_count > 0:
             print("The graph might be disconnected and some vertices are unreachable from S.")
             print("The only path is S itself. Path cover size is 1.")
        return

    # Step 2 & 3: Find the maximum layer size
    layer_sizes = [count for dist, count in sorted(layer_counts.items())]
    
    if len(distances) != V_count:
        print("Warning: Not all vertices are reachable from S.")

    max_layer_size = max(layer_sizes)

    # Output the explanation and result
    layer_sizes_str = ", ".join(map(str, layer_sizes))
    print(f"The sizes of the layers by distance from S are: {layer_sizes_str}")
    print(f"The minimum number of paths required is max({layer_sizes_str}) = {max_layer_size}")


# Example 1:
print("--- Example 1 ---")
# V = {S, 1, 2, 3, 4}, E = {(S, 1),  (S, 2), (1, 3), (2, 3), (2,4)}
# Renaming S to 0 for easier 0-based indexing
V1_count = 5
E1 = [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4)]
S1 = 0
solve_path_cover(V1_count, E1, S1)
# Expected Output: Layers: {0}, {1,2}, {3,4}. Sizes: 1, 2, 2. Max: 2.

print("\n--- Example 2 ---")
# V = {S, 1, ..., 8}, E = {(S, 1), (S, 2), (S, 3), (1, 4), (2, 4), (2, 5), (3, 5), (4, 6), (4, 7), (5, 8)}
# Renaming S to 0
V2_count = 9
E2 = [(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (2, 5), (3, 5), (4, 6), (4, 7), (5, 8)]
S2 = 0
solve_path_cover(V2_count, E2, S2)
# Expected Output: Layers: {0}, {1,2,3}, {4,5}, {6,7,8}. Sizes: 1, 3, 2, 3. Max: 3.
