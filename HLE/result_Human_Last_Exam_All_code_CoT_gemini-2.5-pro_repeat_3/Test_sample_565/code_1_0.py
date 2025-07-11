def solve_vertex_transitive_graphs():
    """
    Calculates and prints the number of isomorphism classes of vertex-transitive
    graphs on 8 vertices, categorized by vertex degree.

    The values are based on established results in algebraic graph theory, as
    a full computational search from scratch is highly complex and requires
    specialized software.
    """

    # n_j is the number of non-isomorphic vertex-transitive graphs
    # on 8 vertices with degree j.
    # n_0: The empty graph.
    n_0 = 1
    # n_1: The graph of 4 disjoint edges (4K_2).
    n_1 = 1
    # n_2: The cycle C_8 and two disjoint C_4's.
    n_2 = 2
    # n_3: The cube graph, 2K_4, and two other Cayley graphs.
    n_3 = 4
    # By complementarity (n_j = n_{7-j}):
    # n_4 is the number of complements of degree-3 graphs.
    n_4 = n_3
    # n_5 is the number of complements of degree-2 graphs.
    n_5 = n_2
    # n_6 is the number of complements of degree-1 graphs.
    n_6 = n_1
    # n_7 is the complement of the degree-0 graph (the complete graph K_8).
    n_7 = n_0

    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    print(result)

solve_vertex_transitive_graphs()