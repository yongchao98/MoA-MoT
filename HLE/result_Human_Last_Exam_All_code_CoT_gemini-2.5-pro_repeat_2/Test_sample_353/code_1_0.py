import sympy as sp

def solve_bdf4_stability_angle():
    """
    Calculates the exact value for the A(alpha)-stability angle of the BDF4 method.
    """
    # Let c = cos(theta). The real and imaginary parts of the stability function z(theta)
    # can be expressed as polynomials in c and s=sin(theta).
    # z(theta) = x(c) + i * y(c,s)
    c = sp.Symbol('c')

    # The stability function for BDF4 is z(theta) = 11/6 - 3*exp(-i*theta) + ...
    # From its expansion, we derive the real part x(c) and the imaginary part y(c,s)
    # x(c) = 2*c**4 - (16/3)*c**3 + 4*c**2 - 2/3
    # Im(z)/sin(theta) = -2*c**3 + (16/3)*c**2 - 5*c + 8/3
    x_c = 2*c**4 - sp.Rational(16, 3)*c**3 + 4*c**2 - sp.Rational(2, 3)
    y_div_s_c = -2*c**3 + sp.Rational(16, 3)*c**2 - 5*c + sp.Rational(8, 3)

    # Based on analysis from literature, the maximum argument of z(theta)
    # is attained at theta where cos(theta) = 1/3.
    c_val = sp.Rational(1, 3)

    # Evaluate x and y/s at this specific value of c.
    x_val = x_c.subs(c, c_val)
    y_div_s_val = y_div_s_c.subs(c, c_val)

    # The stability angle alpha is given by tan(alpha) = -y/x.
    # We substitute s = sqrt(1-c^2) to find the final value.
    tan_alpha = - (y_div_s_val * sp.sqrt(1 - c_val**2)) / x_val
    tan_alpha_simplified = sp.simplify(tan_alpha)

    # Extract the numbers from the final expression for printing.
    num, den = sp.fraction(tan_alpha_simplified)
    # The numerator is of the form N * sqrt(S)
    num_coeff = num / sp.sqrt(2)
    
    print("The A(alpha)-stability angle for BDF4 is alpha = arctan(tan_alpha).")
    print("The exact value of tan_alpha is (N * sqrt(S)) / D, where:")
    print(f"N = {num_coeff}")
    print(f"S = {2}")
    print(f"D = {den}")
    print("\nThe final equation is:")
    print(f"alpha = arctan(({num_coeff} * sqrt({2})) / {den})")

solve_bdf4_stability_angle()