def solve():
    """
    This function determines and prints the subset of integers from {2, 3, 4, 5, 7, 9, 15}
    for which the number of t-omino tilings of an n x n grid is even for any positive integer n.
    """

    # The reasoning for the inclusion of each number is as follows:
    # t is odd (3, 5, 7, 9, 15): A known mathematical theorem confirms the number of tilings is always even
    # for any n x n grid (since n > 1 implies it's not a 1xt rectangle).
    # t=2: For odd n, the number of tilings is 0 (even). For even n=2k, the number of tilings is a multiple of 2^k, hence even.
    # t=4: For odd n, the number is 0 (even). For n = 2 mod 4, the number is even by a symmetry argument.
    # For n = 0 mod 4, the parity is not easily determined. Lacking a proof for all n, we exclude t=4.

    final_subset = [2, 3, 5, 7, 9, 15]

    print("The subset of integers is:")
    print(final_subset)

solve()