import sympy

def solve_bdf4_stability_angle():
    """
    Calculates the exact A(alpha)-stability angle for the BDF4 method.
    The method follows these steps:
    1. Define BDF4 characteristic polynomials rho(xi) and sigma(xi).
    2. Solve the equation rho*sigma' - rho'*sigma = 0 for xi.
    3. Use the complex root xi0 inside the unit circle to find the tangent point z0.
    4. Calculate the angle alpha from z0.
    """

    # Define the symbolic variable
    xi = sympy.Symbol('xi')

    # Step 1: Define BDF4 characteristic polynomials
    # rho(xi) = (25/12)xi^4 - 4xi^3 + 3xi^2 - (4/3)xi + 1/4
    # To avoid floating point issues, we work with rational numbers (fractions).
    rho_poly = sympy.Poly(
        (sympy.S(25)/12)*xi**4 - 4*xi**3 + 3*xi**2 - (sympy.S(4)/3)*xi + sympy.S(1)/4,
        xi
    )
    sigma_poly = sympy.Poly(xi**4, xi)
    
    print("Step 1: The BDF4 characteristic polynomials are:")
    print(f"rho(xi) = {rho_poly.as_expr()}")
    print(f"sigma(xi) = {sigma_poly.as_expr()}")
    print("-" * 30)

    # Step 2: Set up and solve the tangency condition equation
    # The condition is rho(xi)*sigma'(xi) - rho'(xi)*sigma(xi) = 0.
    # For sigma(xi) = xi^k, this simplifies to k*rho(xi) - xi*rho'(xi) = 0. Here k=4.
    tangency_eq = 4 * rho_poly.as_expr() - xi * sympy.diff(rho_poly.as_expr(), xi)
    
    # Simplify the equation. The result is a simple cubic polynomial.
    simplified_eq = sympy.simplify(tangency_eq)
    
    # We solve for the roots of this polynomial. To make it canonical, multiply by -1.
    final_eq = sympy.Poly(-simplified_eq, xi)
    
    print("Step 2: The tangency condition simplifies to the polynomial equation:")
    # Print the equation with integer coefficients
    c = final_eq.all_coeffs()
    print(f"{int(c[0])}*xi^3 - {int(abs(c[1]))}*xi^2 + {int(c[2])}*xi - {int(abs(c[3]))} = 0")
    print("-" * 30)
    
    # Find the roots of the equation
    roots = sympy.solve(final_eq, xi)

    # Select the complex root inside the unit disk with positive imaginary part
    xi0 = None
    for r in roots:
        if sympy.im(r) > 0 and sympy.Abs(r) < 1:
            xi0 = r
            
    print(f"Step 3: Find the relevant root xi0 (complex, |xi0|<1).")
    print(f"The root is xi0 = {xi0}")
    print("-" * 30)

    # Step 3: Compute the point of tangency z0
    print("Step 4: Compute the point z0 = rho(xi0)/sigma(xi0).")
    z0 = sympy.simplify(rho_poly.subs(xi, xi0) / sigma_poly.subs(xi, xi0))
    print(f"The point of tangency is z0 = {z0}")
    print("-" * 30)

    # Step 4: Calculate the angle alpha
    # For A(alpha)-stability, the sector is {z | |arg(-z)| <= alpha}.
    # The tangent line from the origin touches the boundary at z0.
    # So, alpha = |arg(-z0)|.
    # tan(alpha) = |Im(-z0) / Re(-z0)| = | -Im(z0) / -Re(z0) | = Im(z0)/(-Re(z0))
    tan_alpha_val = sympy.im(z0) / (-sympy.re(z0))
    
    print("Step 5: Calculate the angle alpha.")
    print(f"The angle alpha is given by tan(alpha) = Im(z0) / (-Re(z0)).")
    im_z0 = sympy.im(z0)
    re_z0 = sympy.re(z0)
    print(f"tan(alpha) = {im_z0} / ({-re_z0}) = {tan_alpha_val}")
    print("-" * 30)
    
    # Final Result
    final_num = sympy.numer(tan_alpha_val)
    final_den = sympy.denom(tan_alpha_val)
    print("Final Answer:")
    print(f"The exact value for the angle alpha is arctan( Im(z0)/(-Re(z0)) ), which is:")
    print(f"alpha = arctan({final_num}/{final_den})")


solve_bdf4_stability_angle()