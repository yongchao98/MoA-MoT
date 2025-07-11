import sympy

def solve_fixed_point_coupling():
    """
    Calculates and displays the leading order expression for the fixed point coupling u*
    in the phi^4 theory near four dimensions.
    """

    # Step 1: Define the symbols and the one-loop beta function.
    # We work in d = 4 - epsilon dimensions. The dimensionless coupling is 'u'.
    # The beta function β(u) describes the change of the coupling with the energy scale.
    u = sympy.Symbol('u')
    epsilon = sympy.Symbol('epsilon')
    pi = sympy.pi

    # To one-loop order, the beta function for phi^4 theory is:
    # β(u) = -ε*u + (3 / (16*π²)) * u²
    beta_u_expr = -epsilon * u + (3 / (16 * pi**2)) * u**2
    
    print("This script calculates the Wilson-Fisher fixed point for the phi^4 theory.")
    print("The one-loop beta function β(u) in d=4-ε dimensions is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('β(u)'), beta_u_expr), use_unicode=True)
    print("\n" + "="*60 + "\n")

    # Step 2: Find the fixed point by solving β(u*) = 0.
    # Fixed points are the zeros of the beta function. They represent scale-invariant points in the RG flow.
    print("A fixed point u* is found by solving the equation β(u*) = 0:")
    equation = sympy.Eq(beta_u_expr, 0)
    sympy.pprint(equation, use_unicode=True)
    print("\n" + "="*60 + "\n")

    # Step 3: Solve the equation for u.
    # The solutions give the fixed points. We are interested in the non-trivial one.
    solutions = sympy.solve(equation, u)
    
    print(f"The equation has two solutions (fixed points):")
    print(f"1. Gaussian (trivial) fixed point: u* = {solutions[0]}")
    print(f"2. Wilson-Fisher (non-trivial) fixed point: u* = {solutions[1]}")
    print("\n" + "="*60 + "\n")
    
    # Step 4: Present the final expression for the non-trivial fixed point u*.
    # The prompt requests to output each number in the final equation.
    # The non-trivial solution is the Wilson-Fisher fixed point u*.
    u_star = solutions[1]
    
    print("The leading order expression for the Wilson-Fisher fixed point coupling is:")

    # We extract the numerical parts of the expression to display them explicitly.
    coeff_part = u_star.as_coeff_Mul()[0] # Gets the numerical coefficient
    numerator_val = coeff_part.p # Numerator of the coefficient (16)
    denominator_val = coeff_part.q # Denominator of the coefficient (3)

    print(f"u* = ({numerator_val} * π**2 * ε) / {denominator_val}")


if __name__ == "__main__":
    solve_fixed_point_coupling()