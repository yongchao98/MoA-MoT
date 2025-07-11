import collections

def get_connected_components(graph, nodes):
    """Finds the connected components of the subgraph induced by a set of nodes."""
    if not nodes:
        return []
    
    nodes = set(nodes)
    visited = set()
    components = []
    
    for start_node in nodes:
        if start_node not in visited:
            comp = set()
            q = collections.deque([start_node])
            visited.add(start_node)
            comp.add(start_node)
            while q:
                u = q.popleft()
                # Consider only neighbors that are in the allowed node set
                for v in graph.get(u, []):
                    if v in nodes and v not in visited:
                        visited.add(v)
                        comp.add(v)
                        q.append(v)
            components.append(frozenset(comp))
    return components

def partitions_to_sets(p):
    """Converts a list of lists/tuples to a set of frozensets for hashing."""
    return frozenset(frozenset(block) for block in p)

def sets_to_lists(p):
    """Converts a frozenset of frozensets back to a list of lists."""
    return sorted([sorted(list(b)) for b in p])

def meet_G(graph, n, p1, p2):
    """Computes the meet of two connected partitions in P(G,n)."""
    blocks1 = [set(b) for b in p1]
    blocks2 = [set(b) for b in p2]
    
    # 1. Compute meet in the full partition lattice Pi_n
    meet_pi_blocks = []
    for b1 in blocks1:
        for b2 in blocks2:
            intersection = b1.intersection(b2)
            if intersection:
                meet_pi_blocks.append(intersection)

    # 2. Refine the Pi_n meet to get the meet in P(G,n)
    final_meet_blocks = []
    for block in meet_pi_blocks:
        components = get_connected_components(graph, block)
        final_meet_blocks.extend(components)
        
    return partitions_to_sets(final_meet_blocks)

def join_G(graph, n, p1, p2):
    """Computes the join of two connected partitions in P(G,n)."""
    # Join in P(G,n) is the same as in Pi_n.
    # We can find components in an auxiliary graph.
    parent = list(range(n + 1))
    def find_set(v):
        if v == parent[v]:
            return v
        parent[v] = find_set(parent[v])
        return parent[v]
    def unite_sets(a, b):
        a = find_set(a)
        b = find_set(b)
        if a != b:
            parent[b] = a

    all_blocks = list(p1) + list(p2)
    for block in all_blocks:
        if not block: continue
        first_node = list(block)[0]
        for node in block:
            unite_sets(first_node, node)

    components = collections.defaultdict(list)
    for i in range(1, n + 1):
        components[find_set(i)].append(i)
        
    return partitions_to_sets(components.values())

def rank_p(n, p):
    """Computes the rank of a partition."""
    return n - len(p)

def main():
    """
    Demonstrates lattice properties for a specific graph and partitions.
    """
    # Graph G from the analysis: n=6
    # V = {1,2,3,4,5,6}
    # E = {(1,2), (2,3), (1,4), (3,5), (4,6), (5,6)}
    n = 6
    G = {
        1: [2, 4], 2: [1, 3], 3: [2, 5],
        4: [1, 6], 5: [3, 6], 6: [4, 5]
    }

    # Two connected partitions x and y
    x = partitions_to_sets([[1, 2, 3, 5], [4, 6]])
    y = partitions_to_sets([[1, 4, 6, 5], [2, 3]])

    # Compute meet and join in P(G,n)
    meet_xy = meet_G(G, n, x, y)
    join_xy = join_G(G, n, x, y)

    # Compute ranks
    r_x = rank_p(n, x)
    r_y = rank_p(n, y)
    r_meet = rank_p(n, meet_xy)
    r_join = rank_p(n, join_xy)

    print(f"Let G be a graph with V=[6] and E={{1,2},{2,3},{1,4},{3,5},{4,6},{5,6}}.")
    print(f"Let x = {sets_to_lists(x)}")
    print(f"Let y = {sets_to_lists(y)}")
    print("-" * 20)
    print(f"The meet x ∧ y is: {sets_to_lists(meet_xy)}")
    print(f"The join  x ∨ y is: {sets_to_lists(join_xy)}")
    print("-" * 20)
    print("Verifying the semimodular inequality: r(x) + r(y) >= r(x ∧ y) + r(x ∨ y)")
    print(f"r(x) = {n} - {len(x)} = {r_x}")
    print(f"r(y) = {n} - {len(y)} = {r_y}")
    print(f"r(x ∧ y) = {n} - {len(meet_xy)} = {r_meet}")
    print(f"r(x ∨ y) = {n} - {len(join_xy)} = {r_join}")
    print("-" * 20)
    
    lhs = r_x + r_y
    rhs = r_meet + r_join
    
    print(f"r(x) + r(y) = {r_x} + {r_y} = {lhs}")
    print(f"r(x ∧ y) + r(x ∨ y) = {r_meet} + {r_join} = {rhs}")
    print(f"Inequality holds: {lhs} >= {rhs} is {lhs >= rhs}")
    print("\nThis provides computational evidence that the poset is a geometric lattice.")
    print("Combined with the theoretical analysis, the correct statement is B.")


if __name__ == "__main__":
    main()
