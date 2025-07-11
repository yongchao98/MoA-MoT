import sympy as sp

def solve_fixed_point():
    """
    This script calculates the leading order expression for the fixed point coupling u*
    in phi^4 theory near four dimensions (d = 4 - epsilon).
    """

    # 1. Define symbolic variables
    # u_star represents the fixed-point coupling u*
    # epsilon represents the small deviation from 4 dimensions, epsilon = 4 - d
    u_star = sp.Symbol('u^*')
    epsilon = sp.Symbol('epsilon')

    # 2. Define the one-loop beta function for phi^4 theory
    # For a Lagrangian with an interaction term (u/4!)*phi^4, the beta function is:
    # beta(u) = -epsilon*u + (3 * u^2) / (16 * pi^2)
    # The fixed point u* is where beta(u*) = 0.
    beta_function_at_fp = -epsilon * u_star + 3 * u_star**2 / (16 * sp.pi**2)
    
    # 3. Create the equation to solve: beta(u*) = 0
    fixed_point_equation = sp.Eq(beta_function_at_fp, 0)
    
    print("To find the fixed point coupling u^*, we solve the equation beta(u^*) = 0.")
    print("The one-loop beta function equation is:")
    sp.pprint(fixed_point_equation, use_unicode=True)
    print("-" * 50)

    # 4. Solve the equation for u_star
    solutions = sp.solve(fixed_point_equation, u_star)

    # 5. The solutions include the trivial (Gaussian) fixed point u*=0 and the
    #    non-trivial (Wilson-Fisher) fixed point. We are interested in the latter.
    wilson_fisher_fp = None
    for sol in solutions:
        if sol != 0:
            wilson_fisher_fp = sol
            break

    if wilson_fisher_fp is None:
        print("Could not find the non-trivial fixed point.")
        return

    # 6. Present the final result, highlighting each number in the equation.
    print("The non-trivial Wilson-Fisher fixed point is found to be:")
    
    # Deconstruct the expression to present it clearly as requested.
    # The expression is of the form: (A * pi^2 / B) * epsilon
    coefficient = wilson_fisher_fp / epsilon
    numerical_part = coefficient / sp.pi**2
    num, den = numerical_part.as_numer_denom()

    # Create the final formatted equation string
    final_equation_str = f"u^* = ( {num} * \u03c0\u00b2 / {den} ) * \u03b5"

    print("\nThe final equation for the fixed point coupling is:")
    print(final_equation_str)
    
    print("\nOr, using symbolic representation:")
    sp.pprint(sp.Eq(sp.Symbol('u^*'), wilson_fisher_fp), use_unicode=True)


if __name__ == '__main__':
    solve_fixed_point()
