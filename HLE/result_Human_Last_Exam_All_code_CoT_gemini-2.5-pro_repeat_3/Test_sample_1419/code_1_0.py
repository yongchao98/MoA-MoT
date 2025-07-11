import sympy

def find_wilson_fisher_fixed_point():
    """
    Symbolically derives and prints the leading order expression for the
    Wilson-Fisher fixed point coupling in phi^4 theory near d=4 dimensions.
    """
    # Step 1: Define the symbolic variables.
    # u is the dimensionless coupling constant.
    # epsilon represents the deviation from 4 dimensions, where d = 4 - epsilon.
    # pi is the mathematical constant.
    u, epsilon, pi = sympy.symbols('u epsilon pi')
    u_star_sym = sympy.Symbol('u^*')

    # Step 2: Define the one-loop beta function for the single-component phi^4 theory.
    # This is a standard result from Quantum Field Theory in the MS scheme.
    beta_function = -epsilon * u + (3 / (16 * pi**2)) * u**2

    # Step 3: Solve for the fixed points where beta(u) = 0.
    fixed_points = sympy.solve(beta_function, u)

    # The solver returns two solutions: the trivial (Gaussian) fixed point u=0
    # and the non-trivial (Wilson-Fisher) fixed point. We isolate the non-trivial one.
    u_star_expression = None
    for fp in fixed_points:
        if fp != 0:
            u_star_expression = fp
            break

    # Step 4: Programmatically construct and print the final equation.
    # The prompt requires that each number in the final equation be outputted.
    # We will parse the symbolic expression to build a clear output string.
    if u_star_expression:
        # The expression is of the form: (coefficient) * epsilon
        coefficient, _ = u_star_expression.as_coeff_mul(epsilon)
        
        # Get the numerator and denominator of the coefficient
        numer, denom = coefficient.as_numer_denom()
        
        # In the numerator, separate the numerical part from the pi term
        num_coeff, pi_term = numer.as_coeff_mul()
        
        # Get the base and exponent of the pi term (pi^2)
        pi_base, pi_exponent = pi_term.as_base_exp()
        
        print("The leading order expression for the fixed point coupling u* is:")
        # Print the final, formatted equation
        print(f"{u_star_sym} = ({int(num_coeff)} * {pi_base}^{int(pi_exponent)} / {int(denom)}) * {epsilon}")
    else:
        print("Could not determine the non-trivial fixed point.")

if __name__ == '__main__':
    find_wilson_fisher_fixed_point()