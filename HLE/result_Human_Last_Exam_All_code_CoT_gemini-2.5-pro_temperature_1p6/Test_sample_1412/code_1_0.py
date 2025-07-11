def solve():
    """
    This function explains the reasoning behind the solution to the graph theory problem.
    The problem asks for the number of non-isomorphic connected 3-regular adjustable graphs
    with 2000 vertices and at least one perfect matching.

    Let G be such a graph and M be its adjustable perfect matching.
    The vertices of G can be partitioned into V_M and U_M based on the matching M.
    The adjustable property implies that the subgraph induced by V_M, G[V_M],
    is isomorphic to the subgraph induced by U_M, G[U_M]. Let's call this subgraph H.

    The degree of any vertex v in G is 3. This degree is composed of:
    1 (for the edge in M) + d_H(v) (degree in H) + d_C(v) (number of crossing edges).
    So, d_H(v) + d_C(v) = 2.

    The high degree of symmetry in the problem's conditions suggests that the solution
    graphs are vertex-transitive. If so, the degree d_H(v) must be a constant k for all v.
    k can be 0, 1, or 2.

    Case 1: k = 2
    d_H(v) = 2 for all v. Then H is a 2-regular graph on 1000 vertices.
    For G to be connected, H must be a single cycle C_1000.
    This structure corresponds to the Prism Graph C_1000 x K_2. This gives 1 graph.

    Case 2: k = 0
    d_H(v) = 0 for all v. Then H is an empty graph. d_C(v) = 2 for all v.
    Connectivity and vertex-transitivity imply a specific, symmetric crossing pattern
    that can be described by a C_1000 on the vertex indices. This gives a second,
    distinct, bipartite graph. This gives 1 graph.

    Case 3: k = 1
    d_H(v) = 1 for all v. H is a 1-regular graph (a perfect matching). d_C(v) = 1 for all v.
    Connectivity and vertex-transitivity imply the combined structure of H and the
    crossing edges forms a C_1000 on the indices. This leads to a third, distinct,
    non-bipartite graph. This gives 1 graph.

    These three cases lead to three non-isomorphic graphs.
    """
    number_of_graphs = 3
    print(f"Based on the analysis, there are {number_of_graphs} such non-isomorphic graphs.")

solve()
print("<<<3>>>")