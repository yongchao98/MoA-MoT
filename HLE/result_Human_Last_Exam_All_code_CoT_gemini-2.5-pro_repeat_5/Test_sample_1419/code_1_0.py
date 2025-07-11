import sympy

def find_wilson_fisher_fixed_point():
    """
    Calculates and displays the Wilson-Fisher fixed point for phi^4 theory
    in d = 4 - epsilon dimensions using the one-loop beta function.
    """
    # Define the symbolic variables
    # u: the coupling constant
    # epsilon: the deviation from 4 dimensions (d = 4 - epsilon)
    # N: the number of field components (e.g., N=1 for a single scalar field)
    u = sympy.Symbol('u')
    epsilon = sympy.Symbol('epsilon')
    N = sympy.Symbol('N')

    print("The one-loop beta function for the O(N) phi^4 theory is of the form:")
    print("β(u) = -ε*u + C*u^2")
    print("-" * 30)

    # The coefficient C comes from the one-loop diagram calculation
    # C = (N + 8) / (8 * pi^2)
    C_numerator = N + 8
    C_denominator = 8 * sympy.pi**2
    C = C_numerator / C_denominator

    # Define the beta function
    beta_function = -epsilon * u + C * u**2

    print("For the O(N) model, the beta function is:")
    sympy.pprint(beta_function, use_unicode=True)
    print("-" * 30)

    # To find the fixed points u*, we solve β(u*) = 0
    print("Solving β(u) = 0 to find the fixed points...")
    fixed_points = sympy.solve(beta_function, u)

    # The solutions are the trivial Gaussian fixed point (u*=0)
    # and the non-trivial Wilson-Fisher fixed point.
    gaussian_fp = fixed_points[0]
    wilson_fisher_fp = fixed_points[1]

    print(f"The trivial (Gaussian) fixed point is: u* = {gaussian_fp}")
    print("\nThe non-trivial (Wilson-Fisher) fixed point, u*, is the leading order expression:")

    # Print the final equation with its components
    # The expression is u* = (8 * pi**2 * epsilon) / (N + 8)
    num, den = wilson_fisher_fp.as_numer_denom()

    print("\nFinal Equation:")
    print(f"      {sympy.pretty(num)}")
    print(f"u* = {'-' * len(sympy.pretty(den))}")
    print(f"      {sympy.pretty(den)}\n")

    print("Where the numbers and symbols in the equation are:")
    print(f"Numerator constant: 8")
    print(f"Symbol: pi (π)")
    print(f"Symbol: epsilon (ε)")
    print(f"Denominator symbol: N (number of field components)")
    print(f"Denominator constant: 8")

    # For the specific case of a single scalar field (N=1)
    u_star_N1 = wilson_fisher_fp.subs(N, 1)
    print("\nFor the simplest case of a single scalar field (N=1), the expression is:")
    sympy.pprint(u_star_N1, use_unicode=True)


if __name__ == '__main__':
    find_wilson_fisher_fixed_point()