def count_vertex_transitive_graphs():
    """
    Calculates the numbers of isomorphism classes of vertex-transitive graphs
    on 8 vertices, categorized by vertex degree j.
    """

    # Based on mathematical enumeration and known graph classifications.
    
    # n_j is the number of non-isomorphic vertex-transitive graphs of degree j
    # on 8 vertices.

    # n_0: The empty graph on 8 vertices (8 isolated nodes). Unique and VT.
    n_0 = 1

    # n_1: 4 disjoint edges (a perfect matching). Unique up to isomorphism and VT.
    n_1 = 1

    # n_2: C_8 (8-cycle) and 2*C_4 (two disjoint 4-cycles). Both are VT.
    n_2 = 2

    # n_3: Four cubic graphs.
    # 1. 2*K_4 (disconnected)
    # 2. Cube graph Q_3 (connected, bipartite)
    # 3. Cayley graph of Z_8 with generators {1,4,7} (connected, non-bipartite)
    # 4. A Cayley graph of D_4 (connected, non-bipartite)
    n_3 = 4

    # The number of vertex-transitive j-regular graphs on n vertices equals the
    # number of (n-1-j)-regular ones. Here n=8. So, n_j = n_{7-j}.
    n_4 = n_3
    n_5 = n_2
    n_6 = n_1
    n_7 = n_0

    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    
    # Printing the result in the specified format
    print(f"[{result[0]}, {result[1]}, {result[2]}, {result[3]}, {result[4]}, {result[5]}, {result[6]}, {result[7]}]")

count_vertex_transitive_graphs()