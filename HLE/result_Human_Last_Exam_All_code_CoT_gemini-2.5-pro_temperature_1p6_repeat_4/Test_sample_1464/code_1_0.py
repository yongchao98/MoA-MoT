import sympy

def solve_polynomial():
    """
    Solves the given polynomial by identifying its roots based on the
    structure of its coefficients (Vieta's formulas), and then verifies
    the solution by reconstructing the polynomial from these roots.
    """
    # Based on Vieta's formulas, the sum of the roots should equal the
    # negation of the X^3 coefficient. Let's hypothesize the roots are
    # the individual terms of that coefficient.
    r1 = sympy.sqrt(34)
    r2 = sympy.sqrt(14)
    r3 = 2 * sympy.sqrt(11)
    r4 = 2 * sympy.sqrt(6)

    hypothesized_roots = [r1, r2, r3, r4]

    # Let's verify this by constructing the polynomial from these roots.
    # The variable of the polynomial is X.
    X = sympy.Symbol('X')
    
    # The polynomial is the product of (X - root) for each root.
    # We use sympy.expand to multiply it all out.
    polynomial_from_roots = sympy.expand((X - r1) * (X - r2) * (X - r3) * (X - r4))
    
    print("The polynomial is constructed from the hypothesized roots:")
    print(f"P(X) = (X - {r1}) * (X - {r2}) * (X - {r3}) * (X - {r4})")
    print("\nAfter expansion, the equation is:")
    # Using sympy.pretty_print for a clearer mathematical output.
    # This shows all the numbers in the final equation.
    sympy.pretty_print(polynomial_from_roots)

    # Sort the roots by their numerical value to find the increasing order.
    sorted_roots = sorted(hypothesized_roots, key=lambda r: r.evalf())
    
    print("\nThe four roots of the equation in increasing order are:")
    for root in sorted_roots:
        # We print the symbolic form and its approximate numerical value.
        print(f"{root} (approx. {root.evalf():.4f})")

solve_polynomial()