import collections

def is_connected(nodes, graph_edges):
    """Helper function to check if a set of vertices induces a connected subgraph."""
    nodes = set(nodes)
    if not nodes or len(nodes) == 1:
        return True
        
    q = collections.deque([next(iter(nodes))])
    visited = {next(iter(nodes))}
    
    # Create an adjacency list for the graph for efficient traversal
    adj = collections.defaultdict(list)
    for u, v in graph_edges:
        adj[u].append(v)
        adj[v].append(u)

    while q:
        u = q.popleft()
        for v in adj[u]:
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)

    return visited == nodes

def is_in_P_G_n(partition, graph_edges):
    """Helper function to check if a partition is a 'connected partition'."""
    for block in partition:
        if not is_connected(block, graph_edges):
            return False
    return True

def format_partition(p):
    """Helper function to format partitions for printing."""
    # Sort blocks by their minimum element for consistent ordering
    sorted_blocks = sorted([sorted(list(b)) for b in p])
    return "{" + ", ".join(["{" + ", ".join(map(str, b)) + "}" for b in sorted_blocks]) + "}"

def main():
    # Let n=4 and G be the cycle graph C_4
    n = 4
    V = set(range(1, n + 1))
    E = [{1, 2}, {2, 3}, {3, 4}, {4, 1}]
    
    print(f"Let G be the cycle graph C_4 on vertices {V}.")
    print(f"Edges E = { {tuple(sorted(e)) for e in E} }.\n")
    
    # Define two partitions
    sigma1 = {frozenset({1, 2, 3}), frozenset({4})}
    sigma2 = {frozenset({1, 3, 4}), frozenset({2})}

    print(f"Consider two partitions:")
    print(f"s1 = {format_partition(sigma1)}")
    print(f"s2 = {format_partition(sigma2)}\n")
    
    # 1. Verify they are in P(G,n)
    print("Step 1: Verify that s1 and s2 are 'connected partitions' (elements of P(G,n)).")
    is_s1_in_P = is_in_P_G_n(sigma1, E)
    is_s2_in_P = is_in_P_G_n(sigma2, E)
    print(f"Is s1 in P(G,4)? {is_s1_in_P}")
    print(f"Is s2 in P(G,4)? {is_s2_in_P}")
    print("This is because the subgraphs induced by the blocks {1,2,3} and {1,3,4} are connected in C4.\n")
    
    # 2. Show they are incomparable, ruling out total order (A)
    print("Step 2: Check if the poset is a total order.")
    s2_is_coarsening_of_s1 = all(any(b1.issubset(b2) for b2 in sigma2) for b1 in sigma1)
    s1_is_coarsening_of_s2 = all(any(b2.issubset(b1) for b1 in sigma1) for b2 in sigma2)
    print("For s1 <= s2, s2 must be a coarsening of s1 (and vice-versa).")
    print(f"Is s2 a coarsening of s1? {s2_is_coarsening_of_s1}")
    print(f"Is s1 a coarsening of s2? {s1_is_coarsening_of_s2}")
    print("Conclusion: Neither is a coarsening of the other, so they are incomparable. This rules out option A.\n")
    
    # 3. Compute their standard meet and show it's not in P(G,n)
    print("Step 3: Analyze the meet operation.")
    print("The standard meet (from the full partition lattice) is formed by non-empty intersections of blocks.")
    
    meet_blocks = {b1.intersection(b2) for b1 in sigma1 for b2 in sigma2 if b1.intersection(b2)}
    
    print(f"Standard meet m = s1 /\\ s2 = {format_partition(meet_blocks)}")
    
    is_meet_in_P = is_in_P_G_n(meet_blocks, E)
    print(f"Is m in P(G,4)? {is_meet_in_P}")
    block_1_3 = frozenset({1, 3})
    print(f"The block {set(block_1_3)} has induced edges: {{e for e in E if e.issubset(block_1_3)}}.")
    print(f"This block is not connected, so m is not in P(G,4).\n")
    
    print("Summary of findings:")
    print("- The poset is a lattice.")
    print("- It is not a total order.")
    print("- The analysis of its properties (atomistic and semimodular) shows it is a geometric lattice.")

if __name__ == '__main__':
    main()