import sympy
import numpy as np

def count_real_intersections():
    """
    Calculates the number of real intersection points between two quadrics.
    
    The equations are:
    Q1: 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31 = 0
    Q2: 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149 = 0
    """
    
    # Define symbolic variables
    x, y = sympy.symbols('x y')
    
    # Define the two quadric equations
    q1_expr = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    q2_expr = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149
    
    # Treat the expressions as polynomials in y
    p1 = sympy.Poly(q1_expr, y)
    p2 = sympy.Poly(q2_expr, y)
    
    # Compute the resultant with respect to y to eliminate it.
    # The resultant is a polynomial in x whose roots are the x-coordinates
    # of the intersection points.
    resultant_poly = sympy.resultant(p1, p2)
    
    # Get the coefficients of the resultant polynomial
    coeffs = sympy.Poly(resultant_poly, x).all_coeffs()
    
    # Convert sympy coefficients to floating point numbers for numerical root finding
    coeffs_float = [float(c) for c in coeffs]
    
    # Find all roots (real and complex) using numpy
    roots = np.roots(coeffs_float)
    
    # Count the number of roots that are real (i.e., imaginary part is close to zero)
    real_roots_count = 0
    for root in roots:
        if abs(root.imag) < 1e-9:  # A small tolerance for floating point errors
            real_roots_count += 1
            
    # Each real root for x corresponds to a real intersection point.
    print(real_roots_count)

# Execute the function to find and print the number of real intersections
count_real_intersections()