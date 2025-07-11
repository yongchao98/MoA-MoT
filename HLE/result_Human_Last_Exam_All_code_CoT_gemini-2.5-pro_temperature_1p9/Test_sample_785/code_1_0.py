import sys

# Set a higher recursion limit, though not strictly necessary for this iterative solution.
# This is a good practice for potentially deep computational problems.
sys.setrecursionlimit(2000)

def count_orbits():
    """
    This function calculates the number of orbits by solving a change-making problem.
    The problem is equivalent to finding the number of ways to express 1000 as a sum of
    the dimensions of the irreducible representations of the symmetric group S_5.
    """
    
    # The dimensions of the 7 irreducible representations of S_5.
    # There are two irreps of dimension 1 (trivial and sign), two of dimension 4,
    # two of dimension 5, and one of dimension 6.
    irrep_dims = [1, 1, 4, 4, 5, 5, 6]
    
    # The total dimension of the representation space.
    target_dim = 1000
    
    # Print the equation we are solving, as requested.
    # The number of orbits is the number of non-negative integer solutions (n_1, ..., n_7) to:
    equation_str = " + ".join([f"n_{i+1}*{d}" for i, d in enumerate(irrep_dims)]) + f" = {target_dim}"
    print("The problem reduces to finding the number of non-negative integer solutions to the equation:")
    print(equation_str)
    
    # We use dynamic programming to count the number of solutions.
    # dp[i] will store the number of ways to make sum 'i'.
    dp = [0] * (target_dim + 1)
    
    # There is one way to make sum 0 (by choosing no irreps).
    dp[0] = 1
    
    # For each irreducible representation, we update the dp table.
    # This is analogous to having different types of coins.
    for dim in irrep_dims:
        for j in range(dim, target_dim + 1):
            dp[j] += dp[j - j]
            dp[j] += dp[j - dim]
    
    # After iterating through all the 'coins' (irrep dimensions),
    # dp[target_dim] will hold the total number of ways to form the target sum.
    num_orbits = dp[target_dim]
    
    print("\nThe number of distinct irreducible representations for S_5 are:")
    dims_set = sorted(list(set(irrep_dims)))
    for d in dims_set:
      count = irrep_dims.count(d)
      print(f"Dimension {d}: {count} representation(s)")

    print(f"\nThe target dimension is: {target_dim}")
    print("\nThe total number of orbits is:")
    print(num_orbits)
    
    return num_orbits

# Execute the function to find the answer.
final_answer = count_orbits()

# The final answer in the required format
# print(f"<<<{final_answer}>>>")