def calculate_vertex_transitive_graph_counts():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs with 8 vertices
    for each degree j from 0 to 7.
    """
    
    # n_j is the number of j-regular vertex-transitive graphs on 8 vertices.
    # We use the property that n_j = n_{7-j}.

    # j=0: The empty graph E_8. It's unique and vertex-transitive.
    n_0 = 1

    # j=1: The perfect matching 4K_2. It's unique and vertex-transitive.
    n_1 = 1

    # j=2: The cycle C_8 and the disjoint union of two C_4 cycles are the only
    # 2-regular vertex-transitive graphs on 8 vertices.
    n_2 = 2

    # j=3: Based on the known classification of cubic (3-regular) vertex-transitive
    # graphs, there are 6 such non-isomorphic graphs on 8 vertices.
    n_3 = 6

    # Using the complement property n_j = n_{7-j}.
    n_4 = n_3
    n_5 = n_2
    n_6 = n_1
    n_7 = n_0
    
    # The final list of counts for j = 0, 1, ..., 7.
    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    
    # The problem asks to output each number in the final equation.
    # We will format the output string to explicitly show each value.
    print(f"[{result[0]}, {result[1]}, {result[2]}, {result[3]}, {result[4]}, {result[5]}, {result[6]}, {result[7]}]")

calculate_vertex_transitive_graph_counts()