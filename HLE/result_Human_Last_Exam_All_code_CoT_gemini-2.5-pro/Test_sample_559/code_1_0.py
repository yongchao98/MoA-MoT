import sympy as sp

def find_separatrix():
    """
    This function verifies the separatrix for the given system of differential equations
    and prints the final equation.
    """
    # Step 1: Define symbolic variables and functions
    # t is the independent variable (time)
    # u(t) and d(t) are the dependent functions
    t = sp.Symbol('t')
    u = sp.Function('u')(t)
    d = sp.Function('d')(t)

    # Step 2: Define the system of differential equations from the problem
    # u'(t) = (u(t)-1)*u(t)^2
    # d'(t) = 2*d(t)^2 + (-3*u(t)+5*u(t)^2)*d(t) - u(t)*(1-u(t))*u(t)^2
    u_prime_expr = u**2 * (u - 1)
    d_prime_expr = 2*d**2 + (5*u**2 - 3*u)*d - u**3*(1-u)

    # Step 3: Propose the equation for the separatrix
    # Through analysis (as described in the plan), the candidate for the separatrix is d = -u^2.
    separatrix_eq = -u**2

    # Step 4: Verify that the proposed separatrix is a valid trajectory.
    # This is done by showing that if d(t) = -u(t)^2, it satisfies the first ODE.

    # Left-Hand Side (LHS) of the equation for d': d'(t)
    # We compute this by taking the derivative of the separatrix equation w.r.t. t
    # d/dt (separatrix_eq) = d/dt (-u^2) = -2*u * u'(t)
    lhs = sp.diff(separatrix_eq, t)
    # Now, substitute the expression for u'(t) from the system
    lhs = lhs.subs(sp.diff(u, t), u_prime_expr)
    
    # Right-Hand Side (RHS) of the equation for d':
    # We substitute d = -u^2 into the original expression for d'(t).
    rhs = d_prime_expr.subs(d, separatrix_eq)

    # Step 5: Simplify both sides to check for equality.
    lhs_simplified = sp.simplify(lhs)
    rhs_simplified = sp.simplify(rhs)

    print("Verification of the separatrix d = -u^2:")
    print("-" * 40)
    print(f"LHS after substitution (d'/dt): {lhs_simplified}")
    print(f"RHS after substitution: {rhs_simplified}")
    
    if sp.simplify(lhs - rhs) == 0:
        print("\nVerification successful: LHS equals RHS.")
        print("The curve d = -u^2 is a trajectory of the system and a separatrix.")
    else:
        print("\nVerification failed: LHS does not equal RHS.")

    # Step 6: Output the final answer.
    print("\n" + "="*40)
    print("The equation of the separatrix is:")
    # The equation is d = -u^2, which can be written as d + u^2 = 0.
    final_eq_str = "d = -1 * u**2"
    print(final_eq_str)
    
    print("\nNumbers in the final equation (d = a*u^b):")
    print(f"Coefficient 'a': {-1}")
    print(f"Power 'b': {2}")
    print("="*40)

if __name__ == '__main__':
    find_separatrix()