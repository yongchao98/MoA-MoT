import sympy

def display_upper_bound():
    """
    This function prints the derived upper-bound for the given expression.
    """
    # N is a symbolic variable representing the number of nodes.
    N = sympy.symbols('N', positive=True, integer=True)

    # From the derivation, the upper-bound is sqrt(N).
    # This can be written as 1 * sqrt(N).
    factor = 1
    bound_expression = factor * sympy.sqrt(N)
    
    # In the final equation, the number is 1.
    print("The derived upper-bound for ||B * Q_{0, M}||_inf is:")
    print(f"Bound <= {factor} * sqrt(N)")
    print(f"\nThe expression for the upper-bound is {bound_expression}.")
    print(f"The factor of sqrt(N) is {factor}.")

display_upper_bound()