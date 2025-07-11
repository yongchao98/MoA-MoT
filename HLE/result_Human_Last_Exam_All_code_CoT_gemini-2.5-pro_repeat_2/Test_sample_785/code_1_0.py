def solve():
    """
    Calculates the number of orbits for the given group action.
    This corresponds to the number of non-negative integer solutions to
    d1*n1 + d2*n2 + ... = 1000, where di are the dimensions of the
    irreducible representations of the group S_5.
    """
    N = 1000
    # Dimensions of the 7 irreducible representations of S_5
    dims = [1, 1, 4, 4, 5, 5, 6]

    # dp[j] will store the number of ways to make the sum j
    dp = [0] * (N + 1)
    dp[0] = 1

    for d in dims:
        for j in range(d, N + 1):
            dp[j] += dp[j - d]

    # Print the equation being solved
    equation_parts = []
    for i, d in enumerate(dims, 1):
        equation_parts.append(f"{d}*n{i}")
    
    print("The number of orbits is the number of non-negative integer solutions to the equation:")
    print(f"{' + '.join(equation_parts)} = {N}")
    
    result = dp[N]
    print(f"\nThe number of solutions is: {result}")

    # The final answer in the required format
    # print(f"<<<{result}>>>")

solve()
