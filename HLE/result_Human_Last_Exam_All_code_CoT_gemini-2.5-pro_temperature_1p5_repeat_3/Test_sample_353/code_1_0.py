import sympy as sp

def solve_bdf4_stability_angle():
    """
    Calculates the exact value for tan(alpha) for the A(alpha)-stability of BDF4.
    """
    # Define the symbol for c = cos(theta)
    c = sp.Symbol('c')

    # 1. Solve the equation for c = cos(theta_0)
    # This equation arises from the tangency condition for the stability boundary.
    equation = 2*c**2 - 8*c + 5
    sols = sp.solve(equation, c)

    # We need the solution c0 such that |c0| <= 1.
    # The solutions are 2 - sqrt(6)/2 and 2 + sqrt(6)/2.
    # 2 + sqrt(6)/2 is > 1, so we choose the smaller root.
    c0 = sols[0]

    # sin(theta_0), we take the positive root since we can analyze theta in [0, pi]
    s0 = sp.sqrt(1 - c0**2)

    # 2. Define the BDF4 coefficients for rho(xi) = a4*xi^4 + ... + a0
    a4 = sp.Rational(25, 12)
    a3 = -4
    a2 = 3
    a1 = -sp.Rational(4, 3)
    a0 = sp.Rational(1, 4)

    # 3. Express x = Re(z(theta)) and y = Im(z(theta))
    # z(theta) = rho(exp(i*theta)) / exp(i*4*theta)
    # Re(z(theta)) = a4 + a3*cos(theta) + a2*cos(2*theta) + a1*cos(3*theta) + a0*cos(4*theta)
    # Im(z(theta)) = -(a3*sin(theta) + a2*sin(2*theta) + a1*sin(3*theta) + a0*sin(4*theta))
    # Using Chebyshev polynomials T_n(c) = cos(n*theta) and sin(n*theta)/sin(theta) = U_{n-1}(c)
    T = [1, c, 2*c**2 - 1, 4*c**3 - 3*c, 8*c**4 - 8*c**2 + 1]
    U = [0, 1, 2*c, 4*c**2 - 1, 8*c**3 - 4*c]

    x_expr = a4*T[0] + a3*T[1] + a2*T[2] + a1*T[3] + a0*T[4]
    y_expr_div_s = -(a3*U[1] + a2*U[2] + a1*U[3] + a0*U[4])

    # 4. Substitute c = c0 to find x0 and y0
    x0 = x_expr.subs(c, c0)
    y0 = (s0 * y_expr_div_s).subs(c, c0)
    
    # 5. Calculate tan(alpha) = -y0/x0 and simplify
    tan_alpha = -y0 / x0
    tan_alpha_simplified = sp.simplify(tan_alpha)
    
    # 6. Print the result
    # The final angle alpha is given by arctan of this value.
    # The format should be alpha = arctan(...)
    # Let's extract the numbers for the print statement.
    val = tan_alpha_simplified
    
    # The expression is 2*sqrt(2*sqrt(6)-3)
    num1 = val.args[0] # 2
    inner_expr = val.args[1] # sqrt(...)
    num2 = inner_expr.args[0].args[0] # 2
    inner_sqrt_expr = inner_expr.args[0].args[1] # sqrt(6)
    num3 = inner_sqrt_expr.args[0] # 6
    num4 = inner_expr.args[0].args[2] # -3
    
    print("The A(alpha) stability angle for BDF4 is given by alpha = arctan(X), where X is:")
    print(f"{val}")
    print("\nThe equation for alpha is:")
    print(f"alpha = arctan({num1} * sqrt({num2} * sqrt({num3}) - {abs(num4)}))")

solve_bdf4_stability_angle()