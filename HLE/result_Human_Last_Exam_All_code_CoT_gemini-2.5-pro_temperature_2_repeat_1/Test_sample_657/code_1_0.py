import sympy

def solve_controller():
    """
    This function determines and prints the set of all proper stabilizing controllers
    for the given system H1(s) = s / (s^2 - 1).
    """
    # Initialize sympy for symbolic mathematics
    # 's' is the Laplace variable, and 'K' represents the free parameter K(s)
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # The formula for the family of stabilizing controllers H2(s) is derived
    # from the Youla-Kucera Q-parametrization, where Q(s) is replaced by K(s).
    # H2(s) = K(s) * (1 - H1(s)*K(s))^(-1)
    # H2(s) = (s^2 - 1) * K(s) / (s^2 - 1 - s*K(s))

    # Construct the numerator and denominator symbolically
    numerator_expr = (s**2 - 1) * K
    denominator_expr = (s**2 - 1) - s * K

    # Create the controller transfer function
    H2_expr = numerator_expr / denominator_expr

    # Print the final result in a human-readable format.
    # The coefficients (numbers) are explicitly shown in the expression.
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    
    # We use pretty print to format the mathematical expression nicely.
    # It shows all numbers in the equation, for example, the '1' in s^2 - 1.
    print("\n      " + sympy.pretty(numerator_expr, use_unicode=False))
    print("H_2(s) = " + "-" * 30)
    print("      " + sympy.pretty(denominator_expr, use_unicode=False))
    
    print("\nwhere K(s) is any stable and proper rational function.")


if __name__ == '__main__':
    solve_controller()