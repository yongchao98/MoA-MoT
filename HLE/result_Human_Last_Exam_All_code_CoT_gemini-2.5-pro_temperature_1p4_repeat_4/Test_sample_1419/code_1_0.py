import sympy as sp

def find_fixed_point_coupling():
    """
    This function calculates and displays the leading order expression for the
    fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    # Define the symbolic variables for the coupling 'u' and the expansion parameter 'epsilon'
    u, epsilon = sp.symbols('u epsilon')

    # The one-loop beta function for phi^4 theory in d = 4 - epsilon dimensions
    # beta(u) = -epsilon*u + (3*u**2)/(16*pi**2)
    beta_func = -epsilon * u + (3 * u**2) / (16 * sp.pi**2)

    # Fixed points u* are found by solving beta(u) = 0
    fixed_points = sp.solve(beta_func, u)

    # The solutions include the trivial (Gaussian) fixed point u=0.
    # We are interested in the non-trivial (Wilson-Fisher) fixed point.
    u_star_expr = None
    for fp in fixed_points:
        if fp != 0:
            u_star_expr = fp
            break

    print("The fixed point coupling u* is the non-trivial solution to β(u*) = 0.")
    print("The expression for u* is proportional to epsilon: u* = C * epsilon.")
    print("\nBelow are the components of the final equation for u*:")

    # Extract the numerical components of the expression u* = (16 * pi**2 / 3) * epsilon
    numerator_val = 16
    power_of_pi = 2
    denominator_val = 3

    print(f"u* = (numerator * pi ** power_of_pi / denominator) * epsilon")
    print(f"numerator = {numerator_val}")
    print(f"power_of_pi = {power_of_pi}")
    print(f"denominator = {denominator_val}")

    print("\nThus, the final equation for the fixed point coupling u* is:")

    # Use sympy's pretty printing for a clean mathematical representation
    u_star_symbol = sp.Symbol('u^*')
    final_equation = sp.Eq(u_star_symbol, u_star_expr)
    sp.pprint(final_equation, use_unicode=True)

    # For context, also print the approximate numerical value of the coefficient
    coefficient_val = float(u_star_expr / epsilon)
    print(f"\nThe numerical coefficient C = 16 * π² / 3 is approximately {coefficient_val:.4f}.")


find_fixed_point_coupling()