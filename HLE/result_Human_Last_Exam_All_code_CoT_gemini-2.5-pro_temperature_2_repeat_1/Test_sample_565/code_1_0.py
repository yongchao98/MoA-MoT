def solve_vertex_transitive_graphs():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive
    graphs with 8 vertices and vertex degree precisely j for j=0, ..., 7.
    """
    # Based on theoretical enumeration of vertex-transitive graphs on 8 vertices.

    # n_0: The number of 0-regular vertex-transitive graphs on 8 vertices.
    # This is the empty graph E_8.
    n_0 = 1

    # n_1: The number of 1-regular vertex-transitive graphs on 8 vertices.
    # This is the perfect matching graph, 4K_2.
    n_1 = 1

    # n_2: The number of 2-regular vertex-transitive graphs on 8 vertices.
    # These are the cycle C_8 and the disjoint union of two cycles 2C_4.
    n_2 = 2

    # n_3: The number of 3-regular (cubic) vertex-transitive graphs on 8 vertices.
    # These are the Cube graph, the Mobius ladder M_8, the disjoint union 2K_4,
    # and a non-bipartite graph from the quaternion group.
    n_3 = 4

    # For a vertex-transitive graph G of degree j on n vertices, its complement is
    # vertex-transitive of degree n-1-j. This means n_j = n_{7-j}.
    n_4 = n_3
    n_5 = n_2
    n_6 = n_1
    n_7 = n_0

    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    
    # Print the final result in the specified list format.
    print(f"[{result[0]}, {result[1]}, {result[2]}, {result[3]}, {result[4]}, {result[5]}, {result[6]}, {result[7]}]")

solve_vertex_transitive_graphs()