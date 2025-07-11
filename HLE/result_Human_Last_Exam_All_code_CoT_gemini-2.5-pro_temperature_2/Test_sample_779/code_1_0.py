import collections

def is_connected(nodes, graph_edges):
    """Checks if the subgraph induced by a set of vertices is connected using BFS."""
    if not nodes or len(nodes) == 1:
        return True
    
    # Build adjacency list for the subgraph on the given nodes
    adj = collections.defaultdict(list)
    node_set = set(nodes)
    for u, v in graph_edges:
        if u in node_set and v in node_set:
            adj[u].append(v)
            adj[v].append(u)

    # All nodes in the component must be in the adj list if it's connected and not isolated
    if not all(n in adj for n in node_set):
        return False

    # Standard BFS to check connectivity
    start_node = next(iter(nodes))
    queue = collections.deque([start_node])
    visited = {start_node}
    
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                queue.append(v)
                
    return len(visited) == len(nodes)

def is_in_P(partition, n, graph_edges):
    """Checks if a partition is in P(G, n)."""
    # Verify it is a partition of {1, ..., n}
    all_nodes = {node for block in partition for node in block}
    if all_nodes != set(range(1, n + 1)):
        return False
        
    # Verify each block induces a connected subgraph
    for block in partition:
        if not is_connected(block, graph_edges):
            return False
    return True

def is_refinement(p1, p2):
    """Checks if partition p1 is a refinement of partition p2 (p1 <= p2)."""
    for block1 in p1:
        # Each block in p1 must be a subset of some block in p2
        if not any(block1.issubset(block2) for block2 in p2):
            return False
    return True

# --- Main Analysis: A Counterexample to being a Total Order ---
# Let G be a path graph on n=4 vertices: 1-2-3-4
n = 4
graph_edges = [{1, 2}, {2, 3}, {3, 4}]

# Consider two partitions, sigma1 and sigma2
# Use frozensets for blocks to make them hashable and treated as sets
sigma1 = [frozenset({1, 2}), frozenset({3, 4})]
sigma2 = [frozenset({1}), frozenset({2, 3}), frozenset({4})]

print(f"Let G be a path graph on n={n} vertices.")
print(f"Edges: {graph_edges}")
print(f"Let sigma1 = {[set(b) for b in sigma1]}")
print(f"Let sigma2 = {[set(b) for b in sigma2]}")
print("---")

# 1. Verify that both partitions are in P(G, n)
s1_is_valid = is_in_P(sigma1, n, graph_edges)
s2_is_valid = is_in_P(sigma2, n, graph_edges)
print(f"Is sigma1 in P(G,n)? {s1_is_valid}")
print("Check: Block {1,2} is connected by edge {1,2}. Block {3,4} is connected by edge {3,4}.")
print(f"Is sigma2 in P(G,n)? {s2_is_valid}")
print("Check: Block {1} is connected. Block {2,3} is connected by edge {2,3}. Block {4} is connected.")
print("---")

# 2. Check for comparability
s1_le_s2 = is_refinement(sigma1, sigma2)
s2_le_s1 = is_refinement(sigma2, sigma1)

print(f"Is sigma1 a refinement of sigma2 (sigma1 <= sigma2)? {s1_le_s2}")
print(f"Is sigma2 a refinement of sigma1 (sigma2 <= sigma1)? {s2_le_s1}")
print("---")

if not s1_le_s2 and not s2_le_s1:
    print("Conclusion: sigma1 and sigma2 are incomparable.")
    print("Since there exists a pair of incomparable elements, the poset is not a total order.")
    print("This falsifies statement A.")
else:
    print("An error occurred in the logic or the example.")
