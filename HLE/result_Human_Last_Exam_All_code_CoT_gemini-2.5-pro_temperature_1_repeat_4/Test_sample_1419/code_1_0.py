import sympy

def find_fixed_point_expression():
    """
    This function calculates and prints the leading order expression for the
    Wilson-Fisher fixed point coupling (u*) in phi^4 theory in d=4-epsilon
    dimensions.
    """
    # Define the symbolic variables needed for the calculation.
    # u: the dimensionless coupling constant
    # epsilon: the small parameter, epsilon = 4 - d
    # pi: the mathematical constant
    u, epsilon, pi = sympy.symbols('u epsilon pi')

    # The one-loop beta function for phi^4 theory is of the form:
    # beta(u) = -epsilon*u + B*u**2
    # For a single-component scalar field, the constant B = 3 / (8 * pi**2).
    # We define the numerical parts of B to print them later as requested.
    b_numerator = 3
    b_denominator_coeff = 8
    pi_power = 2
    
    B = b_numerator / (b_denominator_coeff * pi**pi_power)

    # Define the beta function expression.
    beta_function = -epsilon * u + B * u**2

    # A fixed point u* is a solution to beta(u*) = 0.
    # We use sympy.solve to find the roots of the beta function.
    fixed_points = sympy.solve(beta_function, u)

    # The solver returns two fixed points:
    # 1. u = 0 (the trivial Gaussian fixed point)
    # 2. u = epsilon / B (the non-trivial Wilson-Fisher fixed point)
    # We are interested in the non-trivial solution.
    u_star = None
    for fp in fixed_points:
        if fp != 0:
            u_star = fp
            break
            
    # Now, we print the final result.
    print("The leading order expression for the fixed point coupling u* in terms of epsilon is:")
    
    # We manually format the string to show the expression clearly.
    # u_star is epsilon/B = epsilon / (3/(8*pi**2)) = (8*pi**2/3)*epsilon
    final_expression_str = f"u* = ({b_denominator_coeff} * pi**{pi_power} / {b_numerator}) * epsilon"
    print(final_expression_str)
    print()

    # As requested, we print each number present in the final equation.
    print("The numbers in the final equation are:")
    print(b_denominator_coeff)
    print(pi_power)
    print(b_numerator)

if __name__ == '__main__':
    find_fixed_point_expression()