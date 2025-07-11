def solve_markov_hitting_time():
    """
    Calculates and prints the formula for the supremum of alpha.
    
    The derivation shows that sup{alpha} is a function of the constant c.
    This function outputs the components of that function as requested.
    """
    
    # The derived formula is sup{alpha} = 1/2 + 2*c.
    # We output the numbers in this final equation.
    
    term1_numerator = 1
    term1_denominator = 2
    term2_coefficient = 2
    
    print("The problem asks for the supremum of alpha such that the alpha-th moment of the hitting time is finite.")
    print("Based on a diffusion approximation of the Markov chain, the result is derived as a function of the constant c.")
    print("\nThe formula for the supremum of alpha is:")
    print(f"sup(alpha) = {term1_numerator}/{term1_denominator} + {term2_coefficient} * c")

# Execute the function to display the answer.
solve_markov_hitting_time()