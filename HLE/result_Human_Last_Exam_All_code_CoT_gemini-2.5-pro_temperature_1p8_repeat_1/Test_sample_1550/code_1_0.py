import sympy

def solve_intersection():
    """
    Calculates the number of real intersection points for the two given quadrics.
    """
    # Define symbols and the two polynomial equations from the problem
    x, y = sympy.symbols('x y')
    
    # Equation 1: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    eq1 = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    
    # Equation 2: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    eq2 = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    # Treat the equations as polynomials in y, with coefficients being functions of x
    poly1_y = sympy.Poly(eq1, y)
    poly2_y = sympy.Poly(eq2, y)

    # Compute the resultant of the two polynomials with respect to y.
    # The resultant is a polynomial in x whose roots correspond to the x-coordinates
    # of the intersection points of the original two curves.
    resultant_poly = sympy.resultant(poly1_y, poly2_y)

    # The final equation we need to solve is resultant_poly = 0.
    # First, let's display the numbers (coefficients) in this final equation.
    coeffs = resultant_poly.all_coeffs()
    degree = resultant_poly.degree()
    
    print("The elimination of y results in a polynomial equation in x.")
    print("The final equation is:")
    
    equation_parts = []
    for i, c in enumerate(coeffs):
        power = degree - i
        # Format the string for each term of the polynomial
        if c != 0:
            term = ""
            if power > 0:
                if c == -1:
                    term += "-"
                elif c != 1:
                    term += f"{c}*"
                term += f"x^{power}" if power > 1 else "x"
            else:
                term += f"{c}"
            equation_parts.append(term)
    
    # Print the equation in a readable format, handling signs properly
    final_equation_str = " + ".join(equation_parts).replace("+ -", "- ")
    print(f"{final_equation_str} = 0\n")

    print("The numbers in this final equation are its coefficients:")
    for i, c in enumerate(coeffs):
        print(f"  Coefficient of x^{degree - i}: {c}")

    # Now, we find the number of real roots of this resultant polynomial.
    # This count gives the number of real x-coordinates for the intersections.
    try:
        real_x_roots = sympy.real_roots(resultant_poly)
        num_real_intersections = len(real_x_roots)
    except Exception as e:
        print(f"\nCould not compute real roots due to an error: {e}")
        # As an alternative, we can try to solve and count real solutions.
        all_roots = sympy.solve(resultant_poly, x)
        num_real_intersections = sum(1 for root in all_roots if root.is_real)

    print(f"\nThe number of real roots for x in this equation is {num_real_intersections}.")
    print(f"Therefore, the number of real intersection points is {num_real_intersections}.")

if __name__ == '__main__':
    solve_intersection()