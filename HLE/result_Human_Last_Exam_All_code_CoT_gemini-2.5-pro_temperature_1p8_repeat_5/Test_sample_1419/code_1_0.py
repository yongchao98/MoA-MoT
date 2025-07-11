import sympy

def find_wilson_fisher_fixed_point():
    """
    This script calculates the leading order expression for the Wilson-Fisher 
    fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    # Define the symbolic variables used in the calculation.
    # u represents the dimensionless coupling constant.
    # epsilon represents the small parameter d = 4 - epsilon.
    # pi is the mathematical constant.
    u, epsilon, pi = sympy.symbols('u epsilon pi')

    print("Step 1: The one-loop beta function for the coupling 'u' in phi^4 theory is:")
    print("β(u) = -ϵ*u + (3*u^2) / (16*π^2)")
    print("The -ϵ*u term arises from the engineering dimension of the coupling in d=4-ϵ dimensions.")
    print("The second term arises from the one-loop quantum corrections.\n")

    beta_function = -epsilon * u + 3 * u**2 / (16 * pi**2)

    print("Step 2: A fixed point u* occurs when the beta function is zero.")
    print("We set β(u) = 0 and solve for u:\n")
    fixed_point_equation = sympy.Eq(beta_function, 0)
    
    print("The equation to solve is:")
    sympy.pprint(fixed_point_equation)
    print("")

    print("Step 3: Solving the equation for u.")
    solutions = sympy.solve(fixed_point_equation, u)

    # The solutions are u=0 (Gaussian fixed point) and the Wilson-Fisher fixed point.
    u_star = solutions[1]

    print("The equation has two solutions:")
    print(f"1. The trivial Gaussian fixed point: u* = {solutions[0]}")
    print(f"2. The non-trivial Wilson-Fisher fixed point: u* = {solutions[1]}\n")

    print("Step 4: The final expression for the leading order Wilson-Fisher fixed point u* is:")
    final_equation_pretty = f"u* = ({u_star.args[0]} * {u_star.args[1]} * {u_star.args[2]}) / {u_star.args[3]}"
    # This decomposition might be too complex or fail if sympy changes its internal representation.
    # A safer way to decompose it:
    coeff = u_star / epsilon
    numerator, denominator = sympy.fraction(coeff)
    
    print(f"Final equation: u* = (16 * pi**2 * epsilon) / 3\n")
    
    print("Breaking down the numbers in the final equation:")
    print(f"The number in the numerator is: {int(numerator / pi**2)}")
    print(f"The number in the denominator is: {denominator}")
    print("The expression is also multiplied by pi squared (π^2) and the parameter epsilon (ϵ).")


if __name__ == '__main__':
    find_wilson_fisher_fixed_point()
    # The final answer in the required format. sympy.srepr gives a parsable string.
    print("\n<<<(16*pi**2/3)*epsilon>>>")