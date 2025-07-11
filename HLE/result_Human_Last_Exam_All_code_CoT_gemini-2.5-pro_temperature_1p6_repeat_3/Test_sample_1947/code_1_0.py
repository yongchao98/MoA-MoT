def solve():
    """
    This function calculates the coefficients for the number of closed tree-like walks of length 6.

    Let the number of closed tree-like walks of length 6 in a simple graph X be written as:
    c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)
    where e is the number of edges, k is the number of K_3 subgraphs, p is the number of P_4 subgraphs.

    A closed walk is tree-like if the set of unique edges traversed forms a tree.
    For a closed walk of length 6, the length is sum(2*k_i) = 6 over the edges of the tree T,
    so sum(k_i)=3, where k_i >= 1.
    The tree T can have 1, 2, or 3 edges.

    - T has 1 edge (P_2, single edge): k_1=3. Contributes to c_1.
      For each edge, 2 walks (starting from either end). c_1 = 2.
    - T has 2 edges (P_3): k_1+k_2=3 -> (1,2) or (2,1). Contributes to c_4.
      For each P_3, there are 8 walks. c_4 = 8.
    - T has 3 edges (P_4 or K_1,3): k_1+k_2+k_3=3 -> k_i=1. Contributes to c_3 and c_5.
      - P_4: 4 walks (starting from any of the 4 vertices). c_3 = 4.
      - K_1,3: 12 walks (6 starting at center, 6 starting at leaves). c_5 = 12.
    - K_3 (triangle) is not a tree. So c_2 = 0.
    """
    
    # The coefficients determined by the combinatorial counting.
    c1 = 2
    c2 = 0
    c3 = 4
    c4 = 8
    c5 = 12
    
    # Print the coefficients in order.
    print(f"{c1}")
    print(f"{c2}")
    print(f"{c3}")
    print(f"{c4}")
    print(f"{c5}")

solve()