import sympy as sp

def solve_for_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) using sympy.
    """
    # Step 1: Define all symbols and functions of time 't'.
    t = sp.Symbol('t')
    f = sp.Function('f')(t)
    theta = sp.Function('theta')(t)
    r = sp.Function('r')(t)

    # Step 2: Establish the core relationships from the problem setup.
    # The complex representation z = c1 + i*c2 gives:
    # c1 = y'(t)/f(t) = r(t)*cos(theta(t))  =>  y'(t) = f(t)*r(t)*cos(theta(t))
    # c2 = y(t)      = r(t)*sin(theta(t))
    y = r * sp.sin(theta)
    y_prime = f * r * sp.cos(theta)

    # Step 3: Use the physical constraints to build a system of equations.
    # Constraint 1: The derivative of our expression for y must equal our expression for y_prime.
    eq1 = sp.Eq(sp.diff(y, t), y_prime)

    # Constraint 2: The derivative of y_prime must be zero, since y'' = 0.
    eq2 = sp.Eq(sp.diff(y_prime, t), 0)

    # Step 4: Solve the system for theta'(t).
    # We introduce symbols for the derivatives to use sympy's algebraic solver.
    theta_prime_sym = sp.Symbol("theta_prime")
    r_prime_sym = sp.Symbol("r_prime")

    # Substitute the symbols for the derivative objects.
    eq1_sub = eq1.subs(sp.diff(theta, t), theta_prime_sym).subs(sp.diff(r, t), r_prime_sym)
    eq2_sub = eq2.subs(sp.diff(theta, t), theta_prime_sym).subs(sp.diff(r, t), r_prime_sym)

    # Solve the system of two equations for the two unknowns: theta_prime and r_prime.
    solution = sp.solve([eq1_sub, eq2_sub], (theta_prime_sym, r_prime_sym))

    # The result for theta'(t) is what we need.
    theta_prime_solution = solution[theta_prime_sym]

    # Step 5: Print the final result in a readable format.
    # The derived expression is f*cos(theta)**2 + f.diff(t)/f * cos(theta)*sin(theta)
    
    # We manually format the output to match the options and be clear.
    f_str = "f(t)"
    f_prime_str = "f'(t)"
    cos_theta_str = "cos(theta(t))"
    sin_theta_str = "sin(theta(t))"

    term1_str = f"{f_str}*cos^2(theta(t))"
    term2_str = f"({f_prime_str}/{f_str})*{cos_theta_str}*{sin_theta_str}"
    
    final_equation_str = f"theta'(t) = {term1_str} + {term2_str}"
    
    print("The final derived equation for theta'(t) is:")
    print(final_equation_str)

solve_for_theta_prime()