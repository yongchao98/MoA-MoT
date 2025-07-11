import collections

def is_connected(nodes, graph_edges):
    """
    Checks if a subgraph induced by a set of nodes is connected using BFS.
    """
    node_set = set(nodes)
    if not node_set:
        return True
    
    q = collections.deque([next(iter(node_set))])
    visited = {next(iter(node_set))}
    
    while q:
        u = q.popleft()
        for v1, v2 in graph_edges:
            # Check neighbors of u
            v = -1
            if v1 == u:
                v = v2
            elif v2 == u:
                v = v1
            
            # If v is a neighbor, in the subgraph, and not visited
            if v != -1 and v in node_set and v not in visited:
                visited.add(v)
                q.append(v)
                
    return visited == node_set

def is_in_P(partition, graph_edges, n):
    """
    Checks if a partition is in P(G, n).
    A partition is in P(G, n) if each of its blocks induces a connected subgraph in G.
    """
    # Check if it's a valid partition of [1, n]
    all_nodes = set()
    for block in partition:
        for node in block:
            all_nodes.add(node)
    if all_nodes != set(range(1, n + 1)):
        return False
        
    # Check connectivity of each block
    for block in partition:
        if not is_connected(block, graph_edges):
            return False
    return True

def format_partition(p):
    """Helper function to format a partition for printing."""
    return "{{{}}}".format(", ".join("{" + ", ".join(map(str, sorted(list(b)))) + "}" for b in sorted(list(p), key=lambda x: min(x))))

def is_refinement(p1, p2):
    """Checks if partition p1 is a refinement of partition p2."""
    for block1 in p1:
        found_super_block = False
        for block2 in p2:
            if block1.issubset(block2):
                found_super_block = True
                break
        if not found_super_block:
            return False
    return True

def main():
    """
    Demonstrates that the poset of connected partitions is not a total order
    using the C4 graph as a counterexample.
    """
    n = 4
    # G is a cycle graph C4 on vertices {1, 2, 3, 4}
    # Edges: {1,2}, {2,3}, {3,4}, {4,1}
    G_edges = {(1, 2), (2, 3), (3, 4), (4, 1)}

    # Two partitions of {1, 2, 3, 4}
    pi1 = {frozenset({1, 2}), frozenset({3, 4})}
    pi2 = {frozenset({1, 4}), frozenset({2, 3})}

    print(f"Let G be a cycle graph on n={n} vertices: {list(G_edges)}")
    print("-" * 40)
    print(f"Consider partition pi1 = {format_partition(pi1)}")
    print(f"Consider partition pi2 = {format_partition(pi2)}")
    print("-" * 40)

    # Verify both partitions are in P(G, n)
    pi1_in_P = is_in_P(pi1, G_edges, n)
    pi2_in_P = is_in_P(pi2, G_edges, n)

    print(f"Is pi1 a connected partition? {pi1_in_P}")
    print(f"Is pi2 a connected partition? {pi2_in_P}")
    print("-" * 40)

    if pi1_in_P and pi2_in_P:
        print("Both partitions are valid elements of the poset P(G, n).")
        
        # Check for comparability
        pi1_le_pi2 = is_refinement(pi1, pi2)
        pi2_le_pi1 = is_refinement(pi2, pi1)
        
        print(f"Is pi1 a refinement of pi2? {pi1_le_pi2}")
        print(f"Is pi2 a refinement of pi1? {pi2_le_pi1}")
        print("-" * 40)
        
        if not pi1_le_pi2 and not pi2_le_pi1:
            print("Conclusion: pi1 and pi2 are incomparable.")
            print("Since there exist incomparable elements, the poset is not a total order.")
        else:
            print("Something is wrong, they should be incomparable.")
    else:
        print("Something is wrong with the setup, both partitions should be valid.")

if __name__ == "__main__":
    main()
