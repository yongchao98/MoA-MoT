import sympy

def solve_bdf4_stability_angle():
    """
    Calculates the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Let c = cos(theta)
    c = sympy.Symbol('c')

    # Real part of z(theta) expressed in terms of c = cos(theta)
    # Re(z) = (1/12) * (24*c**4 - 64*c**3 + 48*c**2 - 8)
    # We derived this from Re(z) = Re(P(w)*conj(w)**4)/12 for w=exp(i*theta)
    x_expr = sympy.Rational(2, 3) * (3*c**4 - 8*c**3 + 6*c**2 - 1)

    # Imaginary part of z(theta) expressed in terms of c and s=sin(theta)
    # From our derivation, y = (1/12)*(96s - 72sc - 64s^3 - 12sc(2c^2-1))
    # where s = sin(theta)
    s = sympy.Symbol('s')
    y_expr_s = sympy.Rational(1, 12) * (96*s - 72*s*c - 64*s**3 - 12*s*c*(2*c**2-1))
    
    # The condition for tangency leads to a polynomial in u=2*cos(theta)
    # 5*u**4 - 32*u**3 + 72*u**2 - 64*u + 16 = 0
    # The roots are u=2 (triple) and u=2/5.
    # The non-trivial root is u=2/5, which means 2*cos(theta) = 2/5
    cos_theta_val = sympy.Rational(1, 5)

    # Substitute c = cos(theta_val) to find the coordinates of the tangency point
    x0 = x_expr.subs(c, cos_theta_val)
    
    # For s, use s**2 = 1 - c**2
    # y expression can be simplified by substituting s**2=1-c**2
    y_expr_c = y_expr_s.subs(s**3, s*(1-c**2))
    # y = (s/12) * (96 - 72c - 64(1-c**2) - 12c(2c**2-1))
    # y = (s/12) * (96 - 72c - 64 + 64c**2 - 24c**3 + 12c)
    # y = (s/12) * (-24c**3 + 64c**2 - 60c + 32)
    y_final_expr = sympy.Rational(1, 3) * s * (-6*c**3 + 16*c**2 - 15*c + 8)

    y0 = y_final_expr.subs(c, cos_theta_val).subs(s, sympy.sqrt(1 - cos_theta_val**2))
    
    # tan(alpha) = -y0/x0
    tan_alpha = -y0 / x0
    
    # Simplify the expression
    tan_alpha_simplified = sympy.simplify(tan_alpha)
    
    num, den = sympy.fraction(tan_alpha_simplified)
    # The expression is of the form A*sqrt(B)/C
    A = num.args[0]
    B = num.args[1].args[0]
    C = den
    
    print("The A(alpha)-stability angle for BDF4 is given by alpha = arctan(tan_alpha)")
    print(f"The exact value of tan_alpha is: {tan_alpha_simplified}")
    print("\nThis can be written as: alpha = arctan(A * sqrt(B) / C)")
    print(f"Where the numbers in the final equation are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")

solve_bdf4_stability_angle()