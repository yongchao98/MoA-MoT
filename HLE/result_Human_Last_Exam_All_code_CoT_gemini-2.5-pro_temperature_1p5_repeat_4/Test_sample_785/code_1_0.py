def count_orbits():
    """
    This function calculates the number of orbits, which has been reduced to a
    combinatorial problem. The problem is equivalent to finding the number of
    non-negative integer solutions to the equation:
    n1*1 + n4*4 + n5*5 + n6*6 = 1000
    where 1, 4, 5, 6 are the dimensions of the irreducible representations of S_5.
    This is a classic "change-making" problem that can be solved with
    dynamic programming.
    """
    
    target_dimension = 1000
    irrep_dims = [1, 4, 5, 6]
    
    print("The final equation is: n_1 * 1 + n_2 * 4 + n_3 * 5 + n_4 * 6 = 1000")
    
    # dp[i] will store the number of ways to form the sum i
    dp = [0] * (target_dimension + 1)
    
    # Base case: There is one way to make sum 0 (by choosing no irreps)
    dp[0] = 1
    
    for dim in irrep_dims:
        for i in range(dim, target_dimension + 1):
            dp[i] += dp[i - dim]
            
    num_orbits = dp[target_dimension]
    print("The number of orbits is:")
    print(num_orbits)

count_orbits()