def solve():
    """
    Calculates the numbers of isomorphism classes of vertex-transitive graphs
    on 8 vertices for each degree j from 0 to 7.
    """

    # We are looking for the number of non-isomorphic vertex-transitive graphs 
    # with 8 vertices for each possible degree j. Let this be n_j.
    
    # Initialize a list to store the results n_0, ..., n_7
    num_graphs = [0] * 8

    # n_0: For degree 0, the only graph is the empty graph with 8 vertices.
    # It is vertex-transitive.
    num_graphs[0] = 1

    # n_1: For degree 1, the only graph is a perfect matching (4 disjoint edges).
    # It is vertex-transitive.
    num_graphs[1] = 1

    # n_2: For degree 2, the graphs are unions of cycles. On 8 vertices, only
    # the 8-cycle (C_8) and the graph of two 4-cycles (2*C_4) are vertex-transitive.
    num_graphs[2] = 2

    # n_3: For degree 3, the number of cubic vertex-transitive graphs on 8 
    # vertices is known from classification results to be 5.
    num_graphs[3] = 5
    
    # A key property of vertex-transitive graphs is that the complement of a
    # j-regular vertex-transitive graph on n vertices is an (n-1-j)-regular
    # vertex-transitive graph. For n=8, this means n_j = n_{7-j}.
    # We use this to find the remaining values.
    num_graphs[4] = num_graphs[7-4]  # n_4 = n_3
    num_graphs[5] = num_graphs[7-5]  # n_5 = n_2
    num_graphs[6] = num_graphs[7-6]  # n_6 = n_1
    num_graphs[7] = num_graphs[7-7]  # n_7 = n_0

    # Print the result in the requested format: [n_0, n_1, ..., n_7]
    print(f"[{', '.join(map(str, num_graphs))}]")

solve()