import sympy

def find_real_intersections():
    """
    Calculates the number of real intersection points for two quadric equations
    using the resultant method.
    """
    x, y = sympy.symbols('x y')

    # Define the two quadric equations from the user input.
    # Equation 1:
    eq1_coeffs = {'A': 164, 'B': -216, 'C': 72, 'D': -16, 'E': 0, 'F': 31}
    p1 = eq1_coeffs['A']*x**2 + eq1_coeffs['B']*x*y + eq1_coeffs['C']*y**2 + \
         eq1_coeffs['D']*x + eq1_coeffs['E']*y + eq1_coeffs['F']
    
    # Equation 2:
    eq2_coeffs = {'A': 864, 'B': -1056, 'C': 324, 'D': -560, 'E': 324, 'F': 149}
    p2 = eq2_coeffs['A']*x**2 + eq2_coeffs['B']*x*y + eq2_coeffs['C']*y**2 + \
         eq2_coeffs['D']*x + eq2_coeffs['E']*y + eq2_coeffs['F']

    print("Finding the number of real intersection points for the quadrics:")
    # Using sympy.Eq to pretty-print the equations with all numbers
    print("Equation 1:", sympy.Eq(p1, 0))
    print("Equation 2:", sympy.Eq(p2, 0))
    print("\n")

    # Treat the polynomials p1 and p2 as polynomials in y
    p1_poly_y = sympy.poly(p1, y)
    p2_poly_y = sympy.poly(p2, y)

    # Compute the resultant with respect to y. This eliminates y and gives a polynomial in x.
    print("Computing the resultant to eliminate y...")
    res_x = sympy.resultant(p1_poly_y, p2_poly_y)
    
    # Find the roots of the resultant polynomial. These are the x-coordinates of the intersections.
    print("Finding the roots of the resultant polynomial R(x) = 0...")
    # Using nroots to find numerical approximations of the roots
    x_roots = sympy.nroots(res_x, n=15)
    
    print("The x-coordinates of the intersection points are:")
    for root in x_roots:
        print(root)
    
    # Check for real roots. A root is real if its imaginary part is close to zero.
    real_x_roots = []
    tolerance = 1e-9
    for r in x_roots:
        if abs(sympy.im(r)) < tolerance:
            real_x_roots.append(sympy.re(r))
            
    print("\nChecking for real solutions...")
    if not real_x_roots:
        print("All x-coordinates are complex.")
        print("Therefore, there are no real (x, y) intersection points.")
        num_real_solutions = 0
    else:
        # If there were real roots, we would proceed to find the corresponding y values.
        # But in this case, all x_roots are complex.
        # This part of the code is for completeness but will not be executed for this problem.
        real_solutions = []
        for x_val in real_x_roots:
            # Substitute x_val into the equations and solve for y
            eq1_at_x = p1.subs(x, x_val)
            eq2_at_x = p2.subs(x, x_val)
            # Find common roots for y
            y_sols = sympy.solve([eq1_at_x, eq2_at_x], y)
            for y_val in y_sols:
                if y_val.is_real:
                    real_solutions.append((x_val, y_val))
        num_real_solutions = len(real_solutions)

    print("\n-------------------------------------------")
    print(f"Number of real intersection points: {num_real_solutions}")
    print("-------------------------------------------")


if __name__ == '__main__':
    find_real_intersections()
