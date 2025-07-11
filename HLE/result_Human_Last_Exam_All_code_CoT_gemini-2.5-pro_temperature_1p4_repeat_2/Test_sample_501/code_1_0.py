import sympy

def polymer_force_law():
    """
    Calculates and displays the force law for a thermally isolated,
    freely jointed polymer chain.
    
    The function uses the symbolic math library `sympy` to represent the physical
    quantities and derive the force by differentiating the energy expression.
    """
    
    # Define the symbolic variables used in the derivation.
    # x: separation of the ends
    # l: length of one segment (strut)
    # n: number of mass points (segments)
    # E0: kinetic energy of the polymer at zero extension (x=0)
    x, l, n, E0 = sympy.symbols('x ell n E(0)')

    # As derived in the plan, the energy E(x) as a function of extension x is:
    energy_expression = E0 * sympy.exp(x**2 / (n**2 * l**2))

    # The force F(x) is the derivative of the energy E(x) with respect to x.
    force_expression = sympy.diff(energy_expression, x)

    # Now, we print the final derived force law clearly.
    # The final equation consists of several parts which we print explicitly.
    print("The derived force law, F(x), between the ends of a thermally isolated polymer chain is:")
    
    # Let's break down the final equation to show each component as requested.
    # F(x) = (coefficient) * (variable part) * (exponential part)
    coefficient = 2 * E0 / (n**2 * l**2)
    variable_part = x
    exponential_part = sympy.exp(x**2 / (n**2 * l**2))

    print(f"\nForce equation broken down into its constituent parts:")
    print(f"F(x) = ( {coefficient} ) * ( {variable_part} ) * ( {exponential_part} )")

    # We can also present the complete equation in a standard mathematical format.
    final_equation = sympy.Eq(sympy.Symbol('F(x)'), force_expression)
    print("\nIn standard notation, the equation is:")
    sympy.pprint(final_equation, use_unicode=True)

# Execute the function to derive and print the force law.
polymer_force_law()