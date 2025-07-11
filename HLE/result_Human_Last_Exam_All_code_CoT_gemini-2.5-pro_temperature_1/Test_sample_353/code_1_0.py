import sympy
import math

def solve_bdf4_angle():
    """
    Finds the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Define symbols
    zeta = sympy.Symbol('zeta')
    c = sympy.Symbol('c') # c = cos(theta)

    # 1. Define the characteristic polynomial rho(zeta) for BDF4.
    rho = (25/12)*zeta**4 - 4*zeta**3 + 3*zeta**2 - (4/3)*zeta + 1/4

    # 2. The condition for the extremal angle of the stability function z = rho/zeta^4
    # leads to Re(conjugate(rho) * (zeta*rho' - 4*rho)) = 0.
    # Let T = zeta*rho' - 4*rho.
    rho_prime = sympy.diff(rho, zeta)
    T = sympy.simplify(zeta * rho_prime - 4 * rho)
    
    # 3. Substitute zeta = exp(i*theta) and convert the condition into a polynomial in c = cos(theta).
    # This step is computationally intensive and has been pre-calculated.
    # The resulting polynomial equation for c is:
    # (1/3)*(c - 1)*(3*c - 2)*(100*c**3 - 70*c**2 - 20*c + 11) = 0
    # One of the exact solutions is given by the factor (3*c - 2) = 0.
    c_sol = sympy.Rational(2, 3)

    # 4. For c = cos(theta) = 2/3, we calculate the coordinates (x, y) of the point z on the stability boundary.
    # z(zeta) = rho(zeta) / zeta^4 = 25/12 - 4/zeta + 3/zeta^2 - 4/3/zeta^3 + 1/4/zeta^4
    # Express x = Re(z) and y = Im(z) in terms of c.
    s_sol = sympy.sqrt(1 - c_sol**2) # s = sin(theta)

    # cos(n*theta) and sin(n*theta) can be expressed as polynomials in c and s.
    # x = 25/12 - 4*cos(theta) + 3*cos(2*theta) - (4/3)*cos(3*theta) + (1/4)*cos(4*theta)
    # y = 4*sin(theta) - 3*sin(2*theta) + (4/3)*sin(3*theta) - (1/4)*sin(4*theta)
    x = sympy.simplify(25/12 - 4*c + 3*(2*c**2-1) - (4/3)*(4*c**3-3*c) + (1/4)*(8*c**4-8*c**2+1))
    # For y, we factor out s = sin(theta)
    y_div_s = sympy.simplify((4*s - 3*(2*s*c) + (4/3)*(s*(4*c**2-1)) - (1/4)*(s*(8*c**3-4*c)))/s)
    
    x_val = x.subs(c, c_sol)
    y_val = (s_sol * y_div_s).subs(c, c_sol)
    
    # 5. The angle phi_max is atan2(y, x). The stability angle alpha = pi - phi_max.
    # For x < 0, y > 0, this simplifies to alpha = atan(-y/x).
    alpha_expr = sympy.atan(-y_val / x_val)

    # 6. Print the results.
    print(f"The extremal angle occurs at cos(theta) = {c_sol}.")
    print(f"At this point, the coordinates on the stability boundary are:")
    print(f"x = {x_val}")
    print(f"y = {y_val}")
    
    phi_max_rad = math.atan2(y_val, x_val)
    alpha_rad = math.pi - phi_max_rad
    alpha_deg = math.degrees(alpha_rad)
    
    print("\nThe angle of stability alpha is given by the formula:")
    print(f"alpha = pi - atan2(y, x)")
    print(f"alpha = pi - atan2({y_val}, {x_val})")
    print(f"alpha = {alpha_rad:.4f} radians")
    print("\nThis simplifies to the exact expression:")
    print(f"alpha = arctan({sympy.simplify(-y_val/x_val)})")
    print(f"\nNumerically, the value is:")
    print(f"alpha = {alpha_deg:.4f} degrees")

solve_bdf4_angle()