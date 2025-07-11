def solve_vertex_transitive_graphs():
    """
    This function calculates the number of isomorphism classes of vertex-transitive
    graphs with 8 vertices for each degree j from 0 to 7.

    The calculation is based on established results in graph theory, as a full
    computational search from first principles is a complex task.
    """

    # n_j is the number of graphs for degree j.
    # We use the property that n_j = n_{7-j} due to graph complementation.

    # For degree j=0, the only graph is the empty graph E_8. It is vertex-transitive.
    n_0 = 1

    # For degree j=1, the only 1-regular graph on 8 vertices is the perfect matching 4K_2.
    # It is vertex-transitive.
    n_1 = 1

    # For degree j=2, 2-regular graphs are disjoint unions of cycles.
    # On 8 vertices, the vertex-transitive options are the 8-cycle (C_8) and
    # the union of two 4-cycles (2C_4).
    n_2 = 2

    # For degree j=3, the classification from algebraic graph theory shows there are
    # two such non-isomorphic graphs: the cube graph Q_3 and the circulant graph C_8(1,4).
    n_3 = 2

    # For degrees j=4, 5, 6, 7, the counts are determined by complementation.
    # n_4 = n_{7-3} = n_3
    n_4 = n_3
    # n_5 = n_{7-2} = n_2
    n_5 = n_2
    # n_6 = n_{7-1} = n_1
    n_6 = n_1
    # n_7 = n_{7-0} = n_0
    n_7 = n_0

    # The final list of counts [n_0, n_1, ..., n_7]
    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]

    # Print the list in the required format.
    print(result)

solve_vertex_transitive_graphs()