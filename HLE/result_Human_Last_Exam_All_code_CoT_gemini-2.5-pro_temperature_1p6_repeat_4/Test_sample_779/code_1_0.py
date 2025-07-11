import collections

def generate_partitions(s):
    """Generates all partitions of a set."""
    if not s:
        yield []
        return
    s = list(s)
    first = s[0]
    rest = s[1:]
    for p in generate_partitions(rest):
        # Option 1: insert first into a new block
        yield [[first]] + p
        # Option 2: insert first into an existing block
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]

def is_connected(nodes, graph_adj):
    """Checks if a subgraph induced by a set of nodes is connected using BFS."""
    if not nodes:
        return True
    nodes = set(nodes)
    q = collections.deque([next(iter(nodes))])
    visited = {next(iter(nodes))}
    
    while q:
        u = q.popleft()
        for v in graph_adj.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    
    return visited == nodes

def is_valid_partition(partition, graph_adj):
    """Checks if a partition is in P(G,n)."""
    for block in partition:
        if not is_connected(block, graph_adj):
            return False
    return True

def is_refinement(p1, p2):
    """Checks if p1 is a refinement of p2 (p1 <= p2)."""
    for block1 in p1:
        found_super_block = False
        for block2 in p2:
            if set(block1).issubset(set(block2)):
                found_super_block = True
                break
        if not found_super_block:
            return False
    return True

def main():
    """
    Finds P(G,n) for a path graph on 4 vertices and shows it's not a total order.
    """
    n = 4
    # G is a path graph 0-1-2-3
    G_adj = {
        0: [1],
        1: [0, 2],
        2: [1, 3],
        3: [2]
    }
    
    vertices = set(range(n))
    all_partitions = list(generate_partitions(vertices))
    
    # P(G,n) is the set of valid partitions
    P_G_n = []
    for p in all_partitions:
        if is_valid_partition(p, G_adj):
            P_G_n.append(p)
            
    print(f"For n={n} and G = Path Graph 0-1-2-3:")
    print(f"The size of P(G,n) is: {len(P_G_n)}")
    
    # Find two incomparable partitions
    for i in range(len(P_G_n)):
        for j in range(i + 1, len(P_G_n)):
            p1 = P_G_n[i]
            p2 = P_G_n[j]
            
            # Sort for consistent output
            p1_sorted = sorted([sorted(b) for b in p1])
            p2_sorted = sorted([sorted(b) for b in p2])

            if not is_refinement(p1, p2) and not is_refinement(p2, p1):
                print("\nFound two incomparable partitions, proving the poset is not a total order:")
                print(f"Partition 1: {p1_sorted}")
                print(f"Partition 2: {p2_sorted}")
                print(f"Is Partition 1 a refinement of Partition 2? {is_refinement(p1, p2)}")
                print(f"Is Partition 2 a refinement of Partition 1? {is_refinement(p2, p1)}")
                return

if __name__ == "__main__":
    main()