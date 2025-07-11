import sympy

def solve_overlap_integral():
    """
    This function calculates and displays the analytical expression for the overlap
    integral of two 2s orbitals in a diatomic molecule like H2+.
    """
    
    # Define the symbolic variables. rho is a dimensionless quantity given by rho = zeta * R,
    # where zeta is the effective nuclear charge and R is the internuclear distance.
    rho = sympy.Symbol('rho')

    # The final analytical expression for the overlap integral S(rho) is of the form:
    # S(rho) = exp(-rho/2) * P(rho)
    # where P(rho) is the following polynomial in rho.
    # The coefficients are derived from the step-by-step integration described in the plan.
    
    c4 = sympy.Rational(1, 240)
    c3 = sympy.Rational(1, 32)
    c2 = sympy.Rational(13, 48)
    c1 = sympy.Rational(5, 4)
    c0 = sympy.Rational(5, 2)

    # Construct the polynomial P(rho)
    polynomial_part = c4 * rho**4 + c3 * rho**3 + c2 * rho**2 + c1 * rho + c0

    # The full expression for S(rho)
    S_rho = sympy.exp(-rho / 2) * polynomial_part

    # Print the explanation and the final formula
    print("The analytical expression for the overlap integral S of two 2s hydrogenic orbitals is a function of rho = zeta * R.")
    print("\nThe full expression is:")
    print(f"S(rho) = exp(-{rho}/2) * P({rho})")
    
    print("\nWhere the polynomial P(rho) is:")
    sympy.pprint(polynomial_part, use_unicode=True)
    
    print("\nBreaking down the final equation into its components:")
    print(f"S(rho) = exp(-{rho}/2) * (c4*{rho}**4 + c3*{rho}**3 + c2*{rho}**2 + c1*{rho} + c0)\n")
    print("The values of the coefficients are:")
    print(f"c4 = 1/240  = {float(c4):.6f}")
    print(f"c3 = 1/32   = {float(c3):.6f}")
    print(f"c2 = 13/48  = {float(c2):.6f}")
    print(f"c1 = 5/4    = {float(c1):.4f}")
    print(f"c0 = 5/2    = {float(c0):.4f}")

    final_expression_str = f"exp(-rho/2) * ( (1/240)*rho**4 + (1/32)*rho**3 + (13/48)*rho**2 + (5/4)*rho + 5/2 )"
    print("\nThus, the final analytical expression is:")
    print(f"S(rho) = {final_expression_str}")

# Execute the function to display the result
solve_overlap_integral()

# The final answer in the requested format
final_answer = "<<<S(rho) = exp(-rho/2) * ( (1/240)*rho**4 + (1/32)*rho**3 + (13/48)*rho**2 + (5/4)*rho + 5/2 )>>>"