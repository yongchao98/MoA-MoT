import sympy

def solve_bdf4_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method
    by finding the point on the stability boundary with the minimum real part.
    """
    # Let u = cos(theta)
    u = sympy.Symbol('u')

    # Chebyshev polynomials T_n(u) represent cos(n*theta)
    # T_0(u) = 1
    # T_1(u) = u
    # T_n(u) = 2*u*T_{n-1}(u) - T_{n-2}(u)
    T = [sympy.Integer(1), u]
    for n in range(2, 5):
        T.append(sympy.simplify(2*u*T[n-1] - T[n-2]))

    # The real part of the BDF4 stability boundary function z(theta) is
    # x(theta) = 25/12 - 4*cos(theta) + 3*cos(2*theta) - (4/3)*cos(3*theta) + (1/4)*cos(4*theta)
    # We express this as a polynomial g(u)
    g = (sympy.Rational(25, 12) * T[0] - 4 * T[1] + 3 * T[2] -
         sympy.Rational(4, 3) * T[3] + sympy.Rational(1, 4) * T[4])
    g = sympy.simplify(g)

    # To find the minimum of g(u) on [-1, 1], we check the endpoints
    # and the critical points (where the derivative is zero).
    g_prime = sympy.diff(g, u)
    critical_points_u = sympy.solve(g_prime, u)

    points_to_check = [sympy.Integer(-1), sympy.Integer(1)]
    for val in critical_points_u:
        # We only care about real roots in the interval [-1, 1]
        if val.is_real and -1 <= val <= 1:
            points_to_check.append(val)

    # Evaluate g(u) at these points to find the minimum value x_m
    min_x = sympy.oo
    min_u = None
    for val in points_to_check:
        current_x = g.subs(u, val)
        if current_x < min_x:
            min_x = current_x
            min_u = val
            
    # For the found u_m = min_u, calculate the corresponding theta and then y_m.
    # min_u is 0, which corresponds to theta = pi/2.
    # y(theta) = 4*sin(theta) - 3*sin(2*theta) + (4/3)*sin(3*theta) - (1/4)*sin(4*theta)
    # At theta = pi/2:
    y_m = (4 * sympy.sin(sympy.pi / 2) - 3 * sympy.sin(sympy.pi) +
           sympy.Rational(4, 3) * sympy.sin(3 * sympy.pi / 2) -
           sympy.Rational(1, 4) * sympy.sin(2 * sympy.pi))
           
    # The point on the boundary is z_m = x_m + i*y_m
    x_m = min_x
    
    # The stability angle alpha is given by arg(-z_m) = arctan(y_m / (-x_m))
    ratio = sympy.simplify(y_m / (-x_m))
    
    # Print the result as an equation showing the numbers involved
    # Final equation is alpha = arctan(y_m / (-x_m))
    print(f"alpha = arctan( ({y_m}) / ({-x_m}) ) = arctan({ratio})")


solve_bdf4_angle()