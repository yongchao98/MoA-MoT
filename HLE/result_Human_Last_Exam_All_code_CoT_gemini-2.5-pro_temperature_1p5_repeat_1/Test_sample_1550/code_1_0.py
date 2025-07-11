import sympy
from math import gcd
from functools import reduce

def count_intersection_points():
    """
    This function calculates the number of real intersection points
    between the two specified quadrics.
    """
    # Define the variables
    x, y = sympy.symbols('x y')

    # Define the two quadric equations as polynomials
    # Eq1: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    # Eq2: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    poly1 = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    poly2 = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    print("The two quadric equations are:")
    print(f"Equation 1: {poly1} = 0")
    print(f"Equation 2: {poly2} = 0")
    print("\n--------------------------------------------------\n")

    # We eliminate 'y' using the resultant method. The resultant is a
    # polynomial in 'x' whose roots correspond to the x-coordinates of
    # the intersection points.
    resultant_poly = sympy.resultant(poly1, poly2, y)

    # The resulting polynomial can have large integer coefficients. We simplify
    # it by creating a polynomial object and dividing by the GCD of its coefficients.
    poly_x = sympy.Poly(resultant_poly, x)
    coeffs = poly_x.all_coeffs()
    
    # In case the coefficients are not integers, we convert them.
    # Here they will be, as they come from integer-coefficient polynomials.
    int_coeffs = [int(c) for c in coeffs]
    common_divisor = reduce(gcd, int_coeffs)
    
    simplified_poly = poly_x.div(common_divisor)[1]

    # Display the final simplified polynomial equation in x.
    # The prompt asks to output each number in the final equation.
    final_coeffs = simplified_poly.all_coeffs()
    c0, c1, c2, c3, c4 = [int(c) for c in final_coeffs]
    
    print("Eliminating 'y' from the system yields a polynomial equation in 'x'.")
    print("The simplified final equation is:")
    print(f"{c0}*x^4 + ({c1})*x^3 + {c2}*x^2 + ({c3})*x + {c4} = 0")
    print("\n--------------------------------------------------\n")

    # The number of real intersection points is the number of real roots of this polynomial.
    # We find the roots of the simplified polynomial.
    roots_x = sympy.solve(simplified_poly, x)

    # Count how many of the roots are real. We check if the imaginary part of the
    # numerical evaluation of each root is zero.
    real_roots_count = 0
    print("The roots of this equation for x are approximately:")
    for root in roots_x:
        root_val = root.evalf()
        print(root_val)
        if sympy.im(root_val) == 0:
            real_roots_count += 1
            
    print("\n--------------------------------------------------\n")
    print(f"The equation has {len(roots_x)} roots in total.")
    print(f"The number of real roots is {real_roots_count}.")
    print("Each real root for 'x' corresponds to a real intersection point.")
    print(f"Therefore, the number of real intersection points is {real_roots_count}.")

if __name__ == '__main__':
    count_intersection_points()