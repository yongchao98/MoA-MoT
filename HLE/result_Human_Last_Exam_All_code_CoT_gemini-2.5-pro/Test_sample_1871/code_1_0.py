import sympy

def solve_derivative():
    """
    This function calculates the partial derivative Dx * rho(alpha, beta)
    by setting up and solving a system of linear equations derived from the
    optimality conditions of the L-infinity distance.
    """
    # Define the symbols for the unknowns:
    # d_rho_d_alpha is the desired derivative, D_x * rho(alpha, beta)
    # d_x0_d_alpha is the derivative of the minimizing x-coordinate w.r.t. alpha
    d_rho_d_alpha = sympy.Symbol('d_rho_d_alpha')
    d_x0_d_alpha = sympy.Symbol('d_x0_d_alpha')

    # Based on our step-by-step derivation, we have the following system of equations:
    # Equation 1: d_rho_d_alpha = d_x0_d_alpha - 1
    # Equation 2: d_rho_d_alpha = -5 * d_x0_d_alpha
    eq1 = sympy.Eq(d_rho_d_alpha, d_x0_d_alpha - 1)
    eq2 = sympy.Eq(d_rho_d_alpha, -5 * d_x0_d_alpha)

    # Solve the system of equations for our two unknowns
    solution = sympy.solve((eq1, eq2), (d_rho_d_alpha, d_x0_d_alpha))

    # The result we are interested in is the value for d_rho_d_alpha
    result = solution[d_rho_d_alpha]

    # The problem asks for the answer as a fraction of two integers.
    # We extract the numerator and denominator from the symbolic fraction.
    if isinstance(result, sympy.Rational):
        numerator = result.p
        denominator = result.q
    else: # Should not happen in this case, but good practice
        frac = sympy.fraction(result)
        numerator = frac[0]
        denominator = frac[1]
        
    print(f"The system of equations to solve is:")
    print(f"1. d(rho)/d(alpha) = d(x0)/d(alpha) - 1")
    print(f"2. d(rho)/d(alpha) = -5 * d(x0)/d(alpha)")
    print("\nSolving this system gives:")
    print(f"The value of D_x * rho(alpha, beta) is {result}.")
    print("\nThe final answer is a fraction of two integers:")
    print(f"The numerator is: {numerator}")
    print(f"The denominator is: {denominator}")
    # Final combined output as requested
    print(f"The final equation is D_x * rho(alpha, beta) = {numerator}/{denominator}")


if __name__ == "__main__":
    solve_derivative()
