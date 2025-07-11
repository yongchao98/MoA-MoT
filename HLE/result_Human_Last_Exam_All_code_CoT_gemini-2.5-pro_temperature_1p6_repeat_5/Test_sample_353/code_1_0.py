import sympy

def solve_bdf4_stability_angle():
    """
    This script calculates the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Define symbolic variables
    xi = sympy.Symbol('xi')
    c = sympy.Symbol('c')
    
    print("Step 1: Define the BDF4 characteristic polynomial rho(xi)")
    # rho(xi) for BDF4
    rho = (25/12)*xi**4 - 4*xi**3 + 3*xi**2 - (4/3)*xi + 1/4
    print("rho(xi) =", rho)
    print("-" * 30)

    print("Step 2: Define the condition for the extremal angle")
    # Condition: Re(xi * rho'(xi) / rho(xi)) = 4 for xi = exp(i*theta)
    # This simplifies to Re(Q(xi) / rho(xi)) = 0, where Q(xi) = xi*rho'(xi) - 4*rho(xi)
    rho_prime = sympy.diff(rho, xi)
    Q = sympy.simplify(xi * rho_prime - 4 * rho)
    print("The condition leads to finding the roots of a polynomial in c = cos(theta).")
    print("Q(xi) = xi*rho'(xi) - 4*rho(xi) =", Q)
    print("-" * 30)
    
    print("Step 3: Formulate and solve the polynomial equation for c = cos(theta)")
    # The condition Re(Q(xi) * conjugate(rho(xi))) = 0 for xi=exp(i*theta)
    # can be transformed into a polynomial in c = cos(theta).
    # This polynomial is 5*c**4 - 16*c**3 + 18*c**2 - 8*c + 1 = 0.
    poly_eqn = 5*c**4 - 16*c**3 + 18*c**2 - 8*c + 1
    print("The polynomial equation for c is:")
    print(f"{poly_eqn} = 0")

    roots = sympy.solve(poly_eqn, c)
    print(f"The roots are: {roots}")
    
    # We are interested in the root that isn't 1 (which corresponds to z=0)
    c0 = [r for r in roots if r != 1][0]
    print(f"The relevant root is c = cos(theta) = {c0}")
    print("-" * 30)

    print("Step 4: Calculate z at the point corresponding to c = cos(theta_0)")
    # s0 = sin(theta_0), choose positive root without loss of generality
    s0 = sympy.sqrt(1 - c0**2)
    print(f"sin(theta_0) = {s0}")

    # Substitute xi = c0 + I*s0 into z(xi) = rho(xi)/xi**4
    sigma = xi**4
    z = rho / sigma
    
    xi0 = c0 + sympy.I * s0
    z0 = z.subs(xi, xi0)
    
    # Separate real and imaginary parts
    re_z0, im_z0 = z0.as_real_imag()
    re_z0 = sympy.simplify(re_z0)
    im_z0 = sympy.simplify(im_z0)
    
    print(f"At this point, z = Re(z) + i*Im(z)")
    print(f"Re(z) = {re_z0}")
    print(f"Im(z) = {im_z0}")
    print("-" * 30)

    print("Step 5: Calculate the stability angle alpha")
    # alpha = arctan(|Im(z)| / |-Re(z)|) since z is in the second quadrant.
    tan_alpha_val = -im_z0 / re_z0
    
    print("alpha = arctan(|Im(z)| / |-Re(z)|)")
    print(f"The argument for arctan is: {tan_alpha_val}")

    # Format the final expression
    num, den = tan_alpha_val.as_numer_denom()
    
    # Check if numerator has sqrt term
    if isinstance(num, sympy.Mul) and any(isinstance(arg, sympy.sqrt) for arg in num.args):
      sqrt_term = [arg for arg in num.args if isinstance(arg, sympy.sqrt)][0]
      coeff = num / sqrt_term
      final_expr_str = f"arctan({coeff}*sqrt({sqrt_term.args[0]})/{den})"
    else:
      final_expr_str = f"arctan({tan_alpha_val})"

    print("\nThe exact value of the angle alpha is:")
    print(final_expr_str)
    
    # Return the value in the specified format
    return f"<<<{final_expr_str}>>>"


final_answer = solve_bdf4_stability_angle()
print(final_answer)
