import networkx as nx

def count_vertex_transitive_graphs():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive graphs
    with 8 vertices and vertex degree precisely j for j=0, ..., 7.
    """
    # Initialize the list of counts [n_0, n_1, ..., n_7]
    n_j = [0] * 8

    # j=0: The only 0-regular graph on 8 vertices is the empty graph. It is vertex-transitive.
    n_j[0] = 1

    # j=1: The only 1-regular graph on 8 vertices is the perfect matching (4 disjoint edges).
    # It is unique up to isomorphism and vertex-transitive.
    n_j[1] = 1

    # j=2: We check all non-isomorphic 2-regular graphs on 8 vertices.
    # The graph6 strings below correspond to C8, 2*C4, and C3+C5.
    g6_2_regular = ["G???_`", "G@?_`?", "GAb`_@"]
    count_2_reg_vt = 0
    for g6 in g6_2_regular:
        g = nx.from_graph6_string(g6)
        if nx.is_vertex_transitive(g):
            count_2_reg_vt += 1
    n_j[2] = count_2_reg_vt

    # j=3: We check all 6 non-isomorphic 3-regular (cubic) graphs on 8 vertices.
    # A cubic graph is vertex-transitive if it's isomorphic to one of the 3 known ones.
    # The known vertex-transitive cubic graphs on 8 vertices are:
    # 1. The Cubical graph
    # 2. The Wagner graph
    # 3. Two disjoint K4's (2*K4)
    vt_cubic_graphs_known = [
        nx.cubical_graph(),
        nx.wagner_graph(),
        nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    ]

    # These are the graph6 strings for all 6 non-isomorphic cubic graphs on 8 vertices.
    g6_cubic = ["G?` D@", "G?a@po", "G?q@`_", "G?cAbo", "G?cAP_", "G@wQ`_"]
    count_3_reg_vt = 0
    for g6 in g6_cubic:
        candidate_graph = nx.from_graph6_string(g6)
        for known_vt_graph in vt_cubic_graphs_known:
            if nx.is_isomorphic(candidate_graph, known_vt_graph):
                count_3_reg_vt += 1
                break
    n_j[3] = count_3_reg_vt

    # For j > 3, we use the complement property: n_j = n_{7-j}.
    n_j[4] = n_j[3]  # n_4 = n_{7-4}
    n_j[5] = n_j[2]  # n_5 = n_{7-5}
    n_j[6] = n_j[1]  # n_6 = n_{7-6}
    n_j[7] = n_j[0]  # n_7 = n_{7-0}

    # Print the final list of counts
    print(n_j)

count_vertex_transitive_graphs()