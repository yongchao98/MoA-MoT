def solve():
    """
    Calculates the number of orbits for the given group action.

    This problem is equivalent to finding the number of ways to express a 1000-dimensional
    representation of the symmetric group S_5 as a direct sum of its irreducible representations.
    This, in turn, is equivalent to finding the number of non-negative integer solutions to
    the Diophantine equation sum(n_i * d_i) = 1000, where d_i are the dimensions of the
    irreducible representations of S_5.
    """

    # The total dimension of the representation space
    target_dimension = 1000

    # The dimensions of the 7 irreducible representations of S_5
    dims = [1, 4, 5, 6, 5, 4, 1]

    # Print the equation we are solving, as requested.
    equation_str = " + ".join([f"n_{i+1}*{d}" for i, d in enumerate(dims)]) + f" = {target_dimension}"
    print("The problem reduces to finding the number of non-negative integer solutions to the equation:")
    print(equation_str)
    
    # dp[k] will store the number of ways to form the sum k
    dp = [0] * (target_dimension + 1)
    dp[0] = 1

    # For each irreducible representation, update the dp table
    for dim in dims:
        for k in range(dim, target_dimension + 1):
            dp[k] += dp[k - dim]

    num_orbits = dp[target_dimension]
    print("\nThe number of orbits is:")
    print(num_orbits)


solve()
<<<58158933350>>>