import sympy as sp

def solve_quadric_intersections():
    """
    Calculates the number of real intersection points for two quadric equations.
    """
    x, y = sp.symbols('x y')

    # The two quadric equations given by the user.
    eq_str1 = "164*x**2 - 216*x*y + 72*y**2 - 16*x + 31"
    eq_str2 = "864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149"

    print(f"Finding the number of real intersection points for the quadrics:")
    print(f"1. {eq_str1.replace('**', '^')} = 0")
    print(f"2. {eq_str2.replace('**', '^')} = 0")
    print("")

    # Convert the string expressions to symbolic expressions
    eq1 = sp.sympify(eq_str1)
    eq2 = sp.sympify(eq_str2)

    # Represent the equations as polynomials in y, with coefficients in x.
    poly1_y = sp.Poly(eq1, y)
    poly2_y = sp.Poly(eq2, y)

    # Compute the resultant of the two polynomials to eliminate y.
    # The resultant is a polynomial in x whose roots are the x-coordinates
    # of the intersection points.
    resultant_poly = sp.resultant(poly1_y, poly2_y)

    # The resultant is a polynomial in x. We create a Poly object for it.
    resultant_poly_x = sp.Poly(resultant_poly, x)
    
    # Find the real roots of the resultant polynomial.
    # The number of distinct real roots corresponds to the number of
    # real x-coordinates of the intersection points.
    real_roots_x = resultant_poly_x.real_roots()

    # The number of real roots gives the number of real intersection points.
    num_intersections = len(real_roots_x)

    print(f"The number of real intersection points is {num_intersections}.")

solve_quadric_intersections()