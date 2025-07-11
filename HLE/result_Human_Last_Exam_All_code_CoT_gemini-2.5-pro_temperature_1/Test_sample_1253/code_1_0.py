import collections

def solve():
    """
    Solves the minimum path cover problem for two example graphs.
    """
    # Example 1
    v1 = ['S', '1', '2', '3', '4']
    e1 = [('S', '1'), ('S', '2'), ('1', '3'), ('2', '3'), ('2', '4')]
    s1 = 'S'
    print("--- Example 1 ---")
    find_min_path_cover(v1, e1, s1)
    print("\n" + "="*20 + "\n")

    # Example 2
    v2 = ['S', '1', '2', '3', '4', '5', '6', '7', '8']
    e2 = [('S', '1'), ('S', '2'), ('S', '3'), ('1', '4'), ('2', '4'), 
          ('2', '5'), ('3', '5'), ('4', '6'), ('4', '7'), ('5', '8')]
    s2 = 'S'
    print("--- Example 2 ---")
    find_min_path_cover(v2, e2, s2)


def find_min_path_cover(vertices, edges, start_node):
    """
    Calculates the size of the smallest set of shortest paths from a source
    that covers all vertices in the graph.

    This is equivalent to finding the size of the largest layer in a BFS traversal.
    """
    if not vertices:
        print("Result = 0")
        return

    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # Layers are stored as a dictionary mapping distance to a list of nodes
    layers = collections.defaultdict(list)
    
    # Standard BFS to find layers
    queue = collections.deque([(start_node, 0)])
    visited = {start_node}
    layers[0].append(start_node)

    while queue:
        current_node, distance = queue.popleft()

        for neighbor in adj[current_node]:
            if neighbor not in visited:
                visited.add(neighbor)
                new_distance = distance + 1
                layers[new_distance].append(neighbor)
                queue.append((neighbor, new_distance))
    
    # Ensure all vertices are covered (graph is connected)
    if len(visited) != len(vertices):
        # This case is not specified in the problem, but a robust solution
        # would handle it. For this problem, we assume a connected graph.
        pass
    
    # Find the maximum size of any layer (excluding the source layer)
    max_size = 0
    layer_sizes = []
    # Sort by distance to print in order
    for dist in sorted(layers.keys()):
        # The source S is always covered by any path, so we can ignore L0
        if dist > 0:
            size = len(layers[dist])
            layer_sizes.append(size)
            if size > max_size:
                max_size = size
    
    # Handling the case of a graph with only the source node
    if not layer_sizes:
        if start_node in vertices:
            max_size = 1 # One path of length 0 to cover S itself
        else:
            max_size = 0

    print(f"Layer sizes (excluding source): {layer_sizes}")
    
    # Building the final equation string as requested
    equation_str = "max(" + ", ".join(map(str, layer_sizes)) + ")" if layer_sizes else "1"
    
    print(f"The size of the smallest set of paths is the maximum of these sizes.")
    print(f"Result = {equation_str} = {max_size}")


# Run the solution
solve()
