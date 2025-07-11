def solve_vertex_transitive_graph_counting():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive
    graphs with 8 vertices and vertex degree precisely j, for j=0 to 7.
    """

    # For a vertex-transitive graph on n=8 vertices, the degree j must be
    # constant for all vertices (the graph is regular).
    #
    # A key property is that if a graph G is vertex-transitive and j-regular,
    # its complement is vertex-transitive and (7-j)-regular. This implies that
    # the number of such graphs for degree j, n_j, is equal to the number
    # for degree 7-j, n_{7-j}.
    #
    # We determine the numbers for j=0, 1, 2, 3 and use symmetry for the rest.

    # n_0: For j=0, the only graph is the empty graph on 8 vertices (E_8).
    # It is vertex-transitive.
    n_0 = 1

    # n_1: For j=1, the only graph is a perfect matching (4K_2).
    # It is vertex-transitive.
    n_1 = 1

    # n_2: For j=2, the graphs are disjoint unions of cycles. On 8 vertices,
    # the vertex-transitive options are the 8-cycle (C_8) and two 4-cycles (2C_4).
    n_2 = 2

    # n_3: For j=3, these are cubic vertex-transitive graphs. From known classifications,
    # there are 4 such graphs on 8 vertices.
    n_3 = 4

    # Using the symmetry property n_j = n_{7-j}:
    n_4 = n_3
    n_5 = n_2
    n_6 = n_1
    n_7 = n_0

    # The final list of numbers.
    result_list = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]

    # Printing the result in the specified format.
    print(f"[{result_list[0]}, {result_list[1]}, {result_list[2]}, {result_list[3]}, {result_list[4]}, {result_list[5]}, {result_list[6]}, {result_list[7]}]")

solve_vertex_transitive_graph_counting()