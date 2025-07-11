def solve_partition_problem():
    """
    Calculates the number of ways to partition a target number into given parts.
    This corresponds to finding the number of non-isomorphic 1000-dimensional
    representations of the group S_5.
    """
    target = 1000
    # Dimensions of the irreducible representations of S_5
    dims = [1, 1, 4, 4, 5, 5, 6]

    # The equation we need to solve is:
    # n1*1 + n2*1 + n3*4 + n4*4 + n5*5 + n6*5 + n7*6 = 1000
    # where ni are non-negative integers.
    
    equation_str = f"{dims[0]}*n1"
    for i in range(1, len(dims)):
        equation_str += f" + {dims[i]}*n{i+1}"
    equation_str += f" = {target}"
    
    print(f"The problem is equivalent to finding the number of non-negative integer solutions to the equation:")
    print(equation_str)
    
    # dp[i] will store the number of ways to make sum i
    dp = [0] * (target + 1)
    dp[0] = 1

    for part in dims:
        for i in range(part, target + 1):
            dp[i] += dp[i - part]

    num_orbits = dp[target]
    print(f"\nThe number of orbits is the number of solutions to this equation, which is:")
    print(num_orbits)

solve_partition_problem()