def solve():
    """
    This function determines the coefficients for the expression of the number of
    closed tree-like walks of length 6 in a simple graph X.

    The expression is of the form:
    c_1*e + c_2*k + c_3*p + c_4*sum(deg(v) choose 2) + c_5*sum(deg(v) choose 3)

    The coefficients are determined by combinatorial counting based on the
    isomorphism class of the tree formed by the edges of the walk.
    """

    # c_1: The tree support is P_2 (an edge).
    # For each edge, there are two walks of length 6 (one starting from each endpoint).
    # u -> v -> u -> v -> u -> v -> u and v -> u -> v -> u -> v -> u -> v.
    c1 = 2

    # c_2: The term corresponds to K_3 (triangles).
    # A tree-like walk requires the subgraph of traversed edges to be a tree.
    # A K_3 is not a tree. Thus, no tree-like walks are counted by this term.
    c2 = 0

    # c_3: The tree support is P_4 (a path on 4 vertices).
    # A detailed combinatorial analysis shows there are 6 such walks for each P_4 subgraph.
    # (1 from each endpoint + 2 from each of the 2 interior vertices)
    c3 = 6

    # c_4: The tree support is P_3 (a path on 3 vertices).
    # A detailed combinatorial analysis shows there are 12 such walks for each P_3 subgraph.
    # (6 from the center vertex + 3 from each of the 2 leaf vertices)
    c4 = 12

    # c_5: The tree support is S_3 (a star graph with 3 edges).
    # A detailed combinatorial analysis shows there are 12 such walks for each S_3 subgraph.
    # (6 from the center vertex + 2 from each of the 3 leaf vertices)
    c5 = 12
    
    # Printing the coefficients in the specified order.
    print(f"{c1},{c2},{c3},{c4},{c5}")

solve()