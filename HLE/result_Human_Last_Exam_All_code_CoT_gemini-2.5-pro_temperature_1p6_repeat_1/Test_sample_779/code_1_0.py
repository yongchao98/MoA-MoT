import collections

def is_connected(nodes, all_edges):
    """
    Checks if the subgraph induced by 'nodes' is connected using BFS.
    The graph is defined by the full list of edges 'all_edges'.
    """
    if not nodes:
        return True
    nodes = list(nodes)
    if len(nodes) == 1:
        return True

    # Build an adjacency list for the subgraph induced by the given nodes.
    adj = collections.defaultdict(list)
    node_set = set(nodes)
    for u, v in all_edges:
        if u in node_set and v in node_set:
            adj[u].append(v)
            adj[v].append(u)
    
    # Check if all nodes are reachable from a single start node.
    start_node = nodes[0]
    q = collections.deque([start_node])
    visited = {start_node}
    
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                q.append(v)
    
    return len(visited) == len(nodes)

def is_in_P_G_n(partition, all_edges):
    """
    Checks if a partition is in P(G,n). This is true if and only if
    every block in the partition induces a connected subgraph in G.
    """
    for block in partition:
        if not is_connected(block, all_edges):
            return False
    return True

def is_refinement(p1, p2):
    """Checks if partition p1 is a refinement of partition p2."""
    for block1 in p1:
        # Each block in p1 must be a subset of some block in p2.
        found_superset = False
        for block2 in p2:
            if block1.issubset(block2):
                found_superset = True
                break
        if not found_superset:
            return False
    return True

def format_partition(p):
    """Helper function for pretty-printing a partition."""
    return str({tuple(sorted(list(b))) for b in p}).replace("'", "")

def main():
    """
    Demonstrates that the poset P(G,n) is not necessarily a total order
    by finding an incomparable pair for a specific graph G.
    """
    n = 4
    # G is the path graph on 4 vertices.
    E = [(1, 2), (2, 3), (3, 4)]
    
    print(f"Analysis of the poset P(G,n) for n={n} and G=Path graph with edges E={E}")

    # Define two partitions. Partitions are represented as frozensets of frozensets.
    sigma1 = frozenset({frozenset({1, 2}), frozenset({3}), frozenset({4})})
    sigma2 = frozenset({frozenset({1}), frozenset({2}), frozenset({3, 4})})

    print(f"\nWe consider two partitions:")
    print(f"  sigma1 = {format_partition(sigma1)}")
    print(f"  sigma2 = {format_partition(sigma2)}")

    # 1. Verify both partitions are in P(G,n).
    print("\nStep 1: Verify that both partitions are in P(G,n).")
    is_s1_valid = is_in_P_G_n(sigma1, E)
    is_s2_valid = is_in_P_G_n(sigma2, E)
    print(f"  Is sigma1 in P(G,n)? {is_s1_valid}. (Block {1,2} is connected by edge (1,2))")
    print(f"  Is sigma2 in P(G,n)? {is_s2_valid}. (Block {3,4} is connected by edge (3,4))")
    
    if not (is_s1_valid and is_s2_valid):
        print("\nError: The chosen partitions are not valid for this example.")
        return

    # 2. Check for the total order property by checking comparability.
    #    The relation sigma_a <= sigma_b holds if sigma_a is a refinement of sigma_b.
    print("\nStep 2: Check for comparability to test the total order property.")
    s1_refines_s2 = is_refinement(sigma1, sigma2)
    print(f"  Is sigma1 <= sigma2 (is sigma1 a refinement of sigma2)? {s1_refines_s2}")
    
    s2_refines_s1 = is_refinement(sigma2, sigma1)
    print(f"  Is sigma2 <= sigma1 (is sigma2 a refinement of sigma1)? {s2_refines_s1}")
    
    # 3. Final conclusion from the script.
    print("\nConclusion from this example:")
    if not s1_refines_s2 and not s2_refines_s1:
        print("  sigma1 and sigma2 are incomparable. Neither is a refinement of the other.")
        print("  Because there exists at least one pair of incomparable elements,")
        print("  the poset is not a total order. This eliminates answer choice A.")
    else:
        print("  This example did not find an incomparable pair.")

if __name__ == '__main__':
    main()