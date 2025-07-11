def solve_and_print():
    """
    Calculates the number of orbits by solving an integer partition problem
    derived from the representation theory of the symmetric group S_5.
    """
    # The dimensions of the 7 irreducible representations of the group S_5.
    irrep_dims = [1, 1, 4, 4, 5, 5, 6]
    
    # The total dimension of the representation space.
    total_dim = 1000

    # The problem is equivalent to finding the number of non-negative integer solutions to
    # the equation: n_1*d_1 + n_2*d_2 + ... + n_7*d_7 = 1000, where d_i are the irrep
    # dimensions and n_i are their multiplicities.

    # We use dynamic programming for this. dp[k] will store the number of ways to form sum k.
    dp = [0] * (total_dim + 1)
    dp[0] = 1

    # For each "coin" (irrep dimension), update the dp table.
    for dim in irrep_dims:
        for j in range(dim, total_dim + 1):
            dp[j] += dp[j - dim]

    # The result is the number of ways to form the total dimension.
    result = dp[total_dim]
    
    # Print the explanation and the final equation as requested.
    print("The number of orbits is equal to the number of non-negative integer solutions to the equation:")
    
    equation_str = (f"n_1 * {irrep_dims[0]} + n_2 * {irrep_dims[1]} + n_3 * {irrep_dims[2]} + "
                    f"n_4 * {irrep_dims[3]} + n_5 * {irrep_dims[4]} + n_6 * {irrep_dims[5]} + "
                    f"n_7 * {irrep_dims[6]} = {total_dim}")
    
    print(equation_str)
    
    print("\nwhere n_i are the multiplicities of the irreducible representations of S_5.")
    
    # Print the final calculated number of orbits.
    print(f"\nThe calculated number of orbits is: {result}")

solve_and_print()