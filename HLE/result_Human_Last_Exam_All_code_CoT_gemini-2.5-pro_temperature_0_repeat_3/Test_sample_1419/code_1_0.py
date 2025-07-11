import sympy as sp

def find_fixed_point_coupling():
    """
    Calculates and displays the leading order expression for the fixed point coupling
    u* in phi^4 theory near four dimensions.
    """
    # Define the symbolic variables for the coupling 'u' and the dimensional parameter 'epsilon'.
    u, epsilon = sp.symbols('u epsilon')
    pi = sp.pi

    print("Step 1: Define the one-loop beta function for phi^4 theory in d = 4 - epsilon dimensions.")
    # The beta function is beta(u) = -epsilon*u + (3 * u^2) / (16 * pi^2)
    # We represent the coefficient as a symbolic fraction for clarity.
    coefficient = sp.Rational(3, 16) * (1 / pi**2)
    beta_function = -epsilon * u + coefficient * u**2
    
    print("The beta function is:")
    # Using sympy's pretty print for a nice mathematical layout.
    sp.pprint(sp.Eq(sp.Symbol('beta(u)'), beta_function), use_unicode=True)
    print("-" * 50)

    print("Step 2: Find the fixed point u* by solving the equation beta(u*) = 0.")
    # A fixed point u* satisfies the condition beta(u*) = 0.
    fixed_point_eq = sp.Eq(beta_function, 0)
    
    print("The equation to solve is:")
    sp.pprint(fixed_point_eq, use_unicode=True)
    print("-" * 50)

    print("Step 3: Solve the equation for u.")
    # Solve the equation for the variable u.
    solutions = sp.solve(fixed_point_eq, u)
    
    # The solutions are the fixed points.
    gaussian_fp = solutions[0]
    wilson_fisher_fp = solutions[1]

    print(f"The equation has two solutions (fixed points):")
    print(f"1. The Gaussian fixed point: u* = {gaussian_fp}")
    print(f"2. The non-trivial Wilson-Fisher fixed point: u* = {wilson_fisher_fp}")
    print("-" * 50)

    print("Step 4: Display the final expression for the non-trivial fixed point coupling u*.")
    # The question asks for the leading order expression for the fixed point coupling u*.
    # This is the non-trivial Wilson-Fisher fixed point.
    
    # Extracting the numbers for the final print statement as requested.
    # u* = (16 * pi**2 / 3) * epsilon
    num_coeff = 16
    den_coeff = 3
    
    print("The leading order expression for the fixed point coupling u* is:")
    print(f"u* = ({num_coeff} * \u03c0\u00b2 / {den_coeff}) * \u03b5")


if __name__ == '__main__':
    find_fixed_point_coupling()