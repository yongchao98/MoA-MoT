def solve():
    """
    This function provides the coefficients for the expression of the number of
    closed tree-like walks of length 6 in a simple graph X.

    The expression is of the form:
    c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)

    The coefficients are derived from a combinatorial analysis of the possible structures
    of such walks. The underlying graph of a tree-like walk must be a tree. For a walk
    of length 6, this tree can have 1, 2, or 3 edges.

    - c1 corresponds to walks on a P_2 (single edge).
    - c2 corresponds to walks on a K_3 (not a tree).
    - c3 corresponds to walks on a P_4 (path of 3 edges).
    - c4 corresponds to walks on a P_3 (path of 2 edges).
    - c5 corresponds to walks on a K_{1,3} (star graph with 3 edges).

    Based on detailed counting for each case, the coefficients are determined.
    """

    # c_1: Coefficient for e (number of edges, P_2)
    c1 = 2

    # c_2: Coefficient for k (number of K_3 subgraphs)
    # The underlying graph of a tree-like walk must be a tree. K_3 is not a tree.
    c2 = 0

    # c_3: Coefficient for p (number of P_4 subgraphs)
    c3 = 6

    # c_4: Coefficient for sum(deg(v) choose 2) (number of P_3 subgraphs)
    c4 = 12

    # c_5: Coefficient for sum(deg(v) choose 3) (number of K_{1,3} subgraphs)
    c5 = 12
    
    # Printing the coefficients in order
    print(c1)
    print(c2)
    print(c3)
    print(c4)
    print(c5)

solve()