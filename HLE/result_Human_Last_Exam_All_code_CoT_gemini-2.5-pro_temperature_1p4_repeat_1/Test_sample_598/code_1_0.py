import networkx as nx

def construct_rook_4x4():
    """
    Constructs the 4x4 Rook's graph, which is an SRG(16,6,2,2).
    This graph is the Cartesian product of two K4 complete graphs.
    """
    k4 = nx.complete_graph(4)
    # The graph is K4 x K4. networkx makes this easy.
    # The nodes will be tuples like (0, 1), (3, 2), etc.
    return nx.cartesian_product(k4, k4)

def construct_shrikhande_graph():
    """
    Constructs the Shrikhande graph, also an SRG(16,6,2,2).
    It can be defined as a Cayley graph on Z_4 x Z_4.
    """
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    
    # The generating set S for the Cayley graph
    # Two vertices v and v+s (mod 4) are connected for s in S.
    S = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    
    for i in range(4):
        for j in range(4):
            for s_i, s_j in S:
                neighbor = ((i + s_i) % 4, (j + s_j) % 4)
                G.add_edge((i, j), neighbor)
    return G

def count_five_cycles(G):
    """
    Counts the number of 5-cycles in a graph G.
    This is a brute-force approach that iterates through all possible paths of length 5.
    Complexity is O(n * d^4) which is feasible for small graphs.
    """
    count = 0
    adj = G.adj
    # We look for a cycle v1-v2-v3-v4-v5-v1
    for v1 in G.nodes():
        for v2 in adj[v1]:
            for v3 in adj[v2]:
                if v3 == v1: continue
                for v4 in adj[v3]:
                    if v4 == v2 or v4 == v1: continue
                    for v5 in adj[v4]:
                        # Path v1-v2-v3-v4-v5 found.
                        # Check if v5 connects back to v1 to close the cycle
                        # and that all vertices in the cycle are distinct.
                        if v5 != v3 and v5 != v2 and v1 in adj[v5]:
                             # At this point v1, v2, v3, v4, v5 are all distinct.
                             # v5 != v1 because it's a neighbor of v1.
                            count += 1
    # Each 5-cycle is counted 10 times: 5 starting nodes and 2 directions.
    return count // 10

if __name__ == "__main__":
    # The parameters for the strongly regular graphs
    n = 16
    d = 6
    l = 2  # lambda
    m = 2  # mu

    # Construct the graphs
    rook_graph = construct_rook_4x4()
    shrikhande_graph = construct_shrikhande_graph()

    # Count the 5-cycles
    c5_rook = count_five_cycles(rook_graph)
    c5_shrikhande = count_five_cycles(shrikhande_graph)

    print(f"Parameters (n, d, lambda, mu): ({n}, {d}, {l}, {m})")
    print("-" * 50)
    print(f"Number of 5-cycles in the Rook's graph (4x4): {c5_rook}")
    print(f"Number of 5-cycles in the Shrikhande graph: {c5_shrikhande}")
    print("-" * 50)
    print("Do they have the same number of 5-cycles?")
    print(f"{c5_rook} == {c5_shrikhande} ? {c5_rook == c5_shrikhande}")
