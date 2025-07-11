def count_orbits():
    """
    This function calculates the number of orbits, which is equivalent to solving
    a partition problem. The problem is to find the number of non-negative integer solutions
    to the equation:
    n_1*d_1 + n_2*d_2 + ... + n_k*d_k = N
    where N=1000 is the dimension of the vector space and d_i are the dimensions
    of the irreducible representations of the group S_5.
    """
    N = 1000
    # The dimensions of the 7 irreducible representations of S_5
    dims = [1, 1, 4, 4, 5, 5, 6]

    # dp[i] will store the number of ways to form the sum i
    dp = [0] * (N + 1)
    dp[0] = 1

    # For each available dimension, update the dp table
    for dim in dims:
        for i in range(dim, N + 1):
            dp[i] += dp[i - dim]

    print("The problem is equivalent to finding the number of non-negative integer solutions to the equation:")
    print("2 * n_1 + 2 * n_4 + 2 * n_5 + 1 * n_6 = 1000, where the coefficients are the number of irreps for each dimension.")
    print("Or more explicitly:")
    print("n_d1_a*1 + n_d1_b*1 + n_d4_a*4 + n_d4_b*4 + n_d5_a*5 + n_d5_b*5 + n_d6_a*6 = 1000")
    print(f"The number of ways to form the sum {N} is: {dp[N]}")
    print(f"<<<{dp[N]}>>>")

count_orbits()