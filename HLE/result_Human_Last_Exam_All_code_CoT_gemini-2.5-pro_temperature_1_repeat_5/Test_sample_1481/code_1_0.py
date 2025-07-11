import sympy

def display_generating_function_asymptotics():
    """
    This function symbolically represents and prints the asymptotic expansion
    of the billiard generating function H(s, s') for small |s' - s|.
    """
    # Define the mathematical symbols
    s, s_prime = sympy.symbols("s s_prime")
    H = sympy.Function('H')(s, s_prime)
    kappa = sympy.Function('κ')(s)

    # Define the terms of the asymptotic expansion
    # Leading term is the absolute arc-length separation
    leading_term = sympy.Abs(s_prime - s)
    
    # First-order correction term due to curvature
    # This term is -(1/24) * κ(s)² * |s' - s|³
    correction_term = (sympy.Integer(1) / sympy.Integer(24)) * (kappa**2) * (sympy.Abs(s_prime - s)**3)

    # Construct the full asymptotic expression
    asymptotic_expression = leading_term - correction_term

    # Create a symbolic equation to represent the approximation
    # H(s, s') ≈ leading_term - correction_term
    asymptotic_equation = sympy.Eq(H, asymptotic_expression, evaluate=False)

    # Print the result in a clear, readable format
    print("The asymptotic expansion of the generating function H(s, s') for a planar Birkhoff billiard is:")
    
    # Using a custom pretty print to ensure the output format
    # The default sympy.pprint can be complex. A formatted string is clearer here.
    lhs_str = f"H(s, s')"
    rhs_str = f"|s' - s| - (1/24)*κ(s)**2*|s' - s|**3"
    
    print(f"\n{lhs_str} ≈ {rhs_str}\n")
    print("Where:")
    print("  s, s' are the arc-length parameters of consecutive boundary collisions.")
    print("  κ(s) is the curvature of the boundary at point s.")
    print("  This approximation holds in the limit |s' - s| → 0.")

if __name__ == '__main__':
    display_generating_function_asymptotics()
