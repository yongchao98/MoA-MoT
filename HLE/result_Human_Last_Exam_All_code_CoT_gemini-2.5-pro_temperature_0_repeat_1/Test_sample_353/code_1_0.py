import sympy

def solve_bdf4_stability_angle():
    """
    Calculates the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Step 1: Define the polynomial equation for x = 2*cos(theta)
    # This equation is derived from the tangency condition for the stability boundary.
    x = sympy.Symbol('x')
    poly_in_x = 5*x**4 - 32*x**3 + 72*x**2 - 64*x + 16

    # Step 2: Solve the polynomial equation for x
    roots_x = sympy.solve(poly_in_x, x)

    # The roots are [2/5, 2]. The root x=2 is trivial (real axis).
    # We need the non-trivial root which defines the tangent.
    x_val = sympy.Rational(2, 5)

    # Step 3: Calculate cos(theta) and sin(theta) from x = 2*cos(theta)
    cos_theta = x_val / 2
    # For the upper-half plane, sin(theta) is positive
    sin_theta = sympy.sqrt(1 - cos_theta**2)

    # Step 4: Calculate the tangency point mu in the complex plane.
    # The stability boundary is given by mu(theta) = P(exp(-i*theta)), where
    # P(w) = 25/12 - 4*w + 3*w^2 - (4/3)*w^3 + (1/4)*w^4.
    w = sympy.Symbol('w')
    P_w = (sympy.Rational(25, 12) - 4*w + 3*w**2 - 
           sympy.Rational(4, 3)*w**3 + sympy.Rational(1, 4)*w**4)
    
    # The point on the unit circle is z_bar = cos(theta) - i*sin(theta)
    z_bar_val = cos_theta - sympy.I * sin_theta
    
    mu_val = sympy.simplify(P_w.subs(w, z_bar_val))

    # Step 5: Extract real and imaginary parts of mu
    x_mu = sympy.re(mu_val)
    y_mu = sympy.im(mu_val)

    # Step 6: Calculate tan(alpha) = -y_mu / x_mu
    tan_alpha = sympy.simplify(-y_mu / x_mu)

    # Step 7: Extract the integer and radical parts for the final expression
    num, den = sympy.fraction(tan_alpha)
    
    # The numerator is of the form C * sqrt(R)
    c, r_sqrt = num.as_coeff_mul()
    rad = r_sqrt[0].args[0] # The number under the square root

    print("The exact value of the angle alpha is given by arctan(C * sqrt(R) / D)")
    print(f"The final equation is alpha = arctan({c} * sqrt({rad}) / {den})")
    print(f"C = {c}")
    print(f"R = {rad}")
    print(f"D = {den}")

solve_bdf4_stability_angle()