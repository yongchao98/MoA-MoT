import sympy

def solve_fixed_point_coupling():
    """
    Calculates and prints the leading order expression for the Wilson-Fisher
    fixed point in phi^4 theory near d=4 dimensions.
    """
    # Define the symbolic variables
    # u_star represents the fixed point coupling constant.
    # epsilon represents the deviation from 4 dimensions, d = 4 - epsilon.
    u_star, epsilon = sympy.symbols('u^* epsilon')
    pi = sympy.pi

    # The one-loop beta function for the dimensionless coupling u in phi^4 theory
    # in d = 4 - epsilon dimensions is given by:
    # β(u) = -ε*u + (3*u^2) / (16*π^2)
    # The first term arises from the canonical dimension of the coupling.
    # The second term is the one-loop quantum correction.
    beta_function_expr = -epsilon * u_star + (3 * u_star**2) / (16 * pi**2)

    print("To find the leading order expression for the fixed point coupling u^* in the ϕ^4 theory,")
    print("we start with the one-loop beta function in d = 4 - ε dimensions:")
    print("β(u) = -ε*u + (3*u^2) / (16*π^2)\n")

    # Fixed points are found by solving β(u^*) = 0.
    # This corresponds to a scale-invariant theory.
    print("A fixed point u^* is a value of the coupling where the theory is scale-invariant,")
    print("which occurs when the beta function is zero. So we solve the equation β(u^*) = 0:")
    print(f"-ε*u^* + (3*(u^*)^2) / (16*π^2) = 0\n")

    # Solve the equation for u_star.
    # sympy.solve finds the values of u_star that make the expression equal to zero.
    solutions = sympy.solve(beta_function_expr, u_star)

    # The equation has two solutions:
    # 1. u^* = 0 (the trivial Gaussian fixed point)
    # 2. The non-trivial Wilson-Fisher fixed point. We select this one.
    wilson_fisher_fp = None
    for sol in solutions:
        if sol != 0:
            wilson_fisher_fp = sol
            break

    print("This equation has a trivial solution u^* = 0 (the Gaussian fixed point) and a non-trivial solution.")
    print("The non-trivial solution is the Wilson-Fisher fixed point, which to leading order in ε is:\n")

    # The solution is u* = (16 * π^2 / 3) * ε.
    # We print each number in the final equation as requested.
    coeff_numerator = 16
    coeff_denominator = 3
    
    print(f"u^* = ({coeff_numerator} * π^2 / {coeff_denominator}) * ε")


if __name__ == '__main__':
    solve_fixed_point_coupling()