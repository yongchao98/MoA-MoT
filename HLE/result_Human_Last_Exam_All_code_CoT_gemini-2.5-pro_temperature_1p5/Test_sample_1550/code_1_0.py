import sympy as sp

def count_real_intersections():
    """
    Calculates the number of real intersection points between two quadrics.
    """
    x, y = sp.symbols('x y', real=True)

    # The first equation: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    eq1_coeffs = [164, -216, 72, -16, 0, 31]
    eq1 = eq1_coeffs[0]*x**2 + eq1_coeffs[1]*x*y + eq1_coeffs[2]*y**2 + eq1_coeffs[3]*x + eq1_coeffs[4]*y + eq1_coeffs[5]
    print("Equation 1: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0")

    # The second equation: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    eq2_coeffs = [864, -1056, 324, -560, 324, 149]
    eq2 = eq2_coeffs[0]*x**2 + eq2_coeffs[1]*x*y + eq2_coeffs[2]*y**2 + eq2_coeffs[3]*x + eq2_coeffs[4]*y + eq2_coeffs[5]
    print("Equation 2: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0")

    # Treat equations as polynomials in y to eliminate y
    p1 = sp.Poly(eq1, y)
    p2 = sp.Poly(eq2, y)

    # Compute the resultant to eliminate y
    res = sp.resultant(p1, p2)
    
    # Find the real roots of the resultant polynomial in x numerically
    # We use a high precision for nroots and a small threshold for checking if imaginary part is zero
    x_roots = [r for r in sp.nroots(res, n=30) if abs(sp.im(r)) < 1e-20]
    
    real_solutions = []
    # Tolerance for checking if two solutions are the same
    tol = 1e-9

    # For each real x, solve for y and check if the solution is consistent
    for x_val in x_roots:
        # Substitute x_val into both equations
        eq1_at_x = eq1.subs(x, x_val)
        eq2_at_x = eq2.subs(x, x_val)
        
        # Find real roots for y from the first equation
        y_sols1 = [r for r in sp.nroots(eq1_at_x, n=30) if abs(sp.im(r)) < 1e-20]
        
        # Check which of these y values also satisfy the second equation
        for y_val in y_sols1:
            # Check if y_val satisfies the second equation within a small tolerance
            if abs(eq2_at_x.subs(y, y_val)) < 1e-9:
                # Check if this solution is new
                is_new = True
                for sol_x, sol_y in real_solutions:
                    if abs(sol_x - x_val) < tol and abs(sol_y - y_val) < tol:
                        is_new = False
                        break
                if is_new:
                    real_solutions.append((x_val, y_val))

    num_real_intersections = len(real_solutions)
    
    print(f"\nThe number of real intersection points is: {num_real_intersections}")

count_real_intersections()