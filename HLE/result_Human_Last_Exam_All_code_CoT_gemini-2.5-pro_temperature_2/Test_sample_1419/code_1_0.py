import sympy as sp

def solve_fixed_point_coupling():
    """
    This function calculates and displays the leading order expression for the
    Wilson-Fisher fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    # Define the symbols for the calculation
    u_star = sp.Symbol('u^*')
    epsilon = sp.Symbol('epsilon')
    pi = sp.pi

    # In phi^4 theory in d = 4 - epsilon dimensions, the one-loop beta function
    # for the dimensionless coupling 'u' is beta(u) = -epsilon*u + B*u^2.
    # The constant B for the single-component theory is 3 / (8*pi^2).
    B = 3 / (8 * pi**2)

    # The fixed point u* is found by solving beta(u*) = 0.
    # We construct the equation: -epsilon*u* + B*(u*)^2 = 0
    fixed_point_equation = sp.Eq(-epsilon * u_star + B * u_star**2, 0)

    # Solve the equation for u*. The result will contain the trivial
    # Gaussian fixed point (u* = 0) and the non-trivial Wilson-Fisher fixed point.
    solutions = sp.solve(fixed_point_equation, u_star)

    # We are interested in the non-trivial solution (where u* is not 0).
    wilson_fisher_fp = None
    for sol in solutions:
        if sol != 0:
            wilson_fisher_fp = sol
            break
            
    if wilson_fisher_fp is None:
        print("Could not find the non-trivial fixed point.")
        return

    # Now, we print the results in a step-by-step fashion.
    print("The beta function at the fixed point u* is set to zero:")
    print(f"β(u*) = -{epsilon}*{u_star} + (3 / (8*{sp.Symbol('pi')}^2)) * {u_star}^2 = 0")
    print("\nSolving for the non-trivial Wilson-Fisher fixed point yields:")
    print(f"{u_star} = {wilson_fisher_fp}")
    print("\nThis is the leading order expression for the fixed point coupling.")
    
    # As requested, outputting each number in the final equation.
    # The solution is of the form: u* = (numerator / denominator) * epsilon
    num_val = 8
    pi_sym = "pi"
    denom_val = 3
    eps_sym = "epsilon"
    
    print("\nFinal equation breakdown:")
    print(f"u* = ({num_val} * {pi_sym}^2 / {denom_val}) * {eps_sym}")

if __name__ == "__main__":
    solve_fixed_point_coupling()