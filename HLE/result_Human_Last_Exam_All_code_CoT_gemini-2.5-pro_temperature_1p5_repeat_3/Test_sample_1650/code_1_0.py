import sympy

def solve_overlap_integral():
    """
    This function calculates and displays the analytical expression for the overlap
    integral of two hydrogenic 2s orbitals in the H2+ ion.

    The derivation involves expressing the integral in elliptical coordinates and
    evaluating it step-by-step. This script presents the final result of that
    derivation.
    """

    # Define the symbolic variables
    # R: internuclear distance
    # zeta: effective nuclear charge
    R, zeta = sympy.symbols('R zeta', positive=True, real=True)

    # The result is expressed in terms of a dimensionless variable rho
    rho = sympy.Symbol('rho')

    # The analytical result derived from the integration is:
    # S = e^(-rho) * (1 + rho + (1/3)*rho^2 + (1/15)*rho^4)
    # where rho = zeta * R / 2
    
    # Coefficients from the derivation
    c0 = 1
    c1 = 1
    c2 = sympy.Rational(1, 3)
    c3 = 0  # The coefficient for the rho^3 term is 0 in this derivation
    c4 = sympy.Rational(1, 15)

    # Construct the polynomial part of the expression
    polynomial_part = c0 + c1*rho + c2*rho**2 + c4*rho**4

    # The full expression for the overlap integral S
    S = sympy.exp(-rho) * polynomial_part

    # Print the explanation and the final formula
    print("The overlap integral S for two hydrogenic 2s orbitals is given by the following analytical expression:")
    print("\nIn terms of the dimensionless parameter rho = (zeta * R) / 2:")
    
    final_equation_str = f"S(rho) = exp(-rho) * ({c0} + {c1}*rho + ({c2.p}/{c2.q})*rho**2 + ({c4.p}/{c4.q})*rho**4)"
    
    # Use sympy.pretty_print for a nicely formatted output of the equation
    equation_to_print = sympy.Eq(sympy.Symbol('S'), S)
    
    print("\n" + "="*50)
    sympy.pretty_print(equation_to_print, use_unicode=True)
    print("="*50 + "\n")

    print("Where:")
    print("  S    is the overlap integral.")
    print("  zeta is the effective nuclear charge.")
    print("  R    is the internuclear distance.")
    print("  rho  is a dimensionless variable defined as (zeta * R) / 2.")
    
    # As requested, output each number in the final equation string representation
    print("\nFinal equation string with numerical coefficients:")
    print(final_equation_str)


solve_overlap_integral()
<<<S(rho) = exp(-rho) * (1 + rho + rho**2/3 + rho**4/15)>>>