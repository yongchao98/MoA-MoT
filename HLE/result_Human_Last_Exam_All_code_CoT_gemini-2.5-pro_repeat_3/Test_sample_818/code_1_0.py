def solve_and_print():
    """
    This function calculates M(0), M(3), and M(5) based on graph theory principles
    and then prints the result in the specified format.
    """

    # M(0): Smallest m for a cubic graph G where N(G) is a multiple of 0 (i.e., N(G) = 0).
    # N(G) = 0 if and only if G does not have a C4-factor.
    # This requires m to NOT be a multiple of 4.
    # The smallest even m (for cubic graphs) not divisible by 4 is 6.
    m_0 = 6

    # M(3): Smallest m where N(G) is a multiple of 3.
    # The smallest cubic graph has m=4 (K4).
    # For K4, the number of slices is N(K4) = 9.
    # 9 is a multiple of 3. Since m=4 is the minimum possible, M(3)=4.
    n_k4 = 9
    if n_k4 % 3 == 0:
        m_3 = 4
    else:
        # This case won't be reached if our N(K4) is correct.
        # If N(K4) was not a multiple of 3, M(3) would be 6 (since N(G)=0 for m=6).
        m_3 = 6


    # M(5): Smallest m where N(G) is a multiple of 5.
    # For m=4, N(K4) = 9, which is not a multiple of 5.
    # For m=6, any cubic graph G has N(G) = 0.
    # 0 is a multiple of 5.
    # So the smallest m is 6.
    if n_k4 % 5 != 0:
        m_5 = 6
    else:
        # If N(K4) was a multiple of 5, M(5) would be 4.
        m_5 = 4

    # Print the final answer in the required format.
    print(f"{m_0},{m_3},{m_5}")

solve_and_print()