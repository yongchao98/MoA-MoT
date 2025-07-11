import sympy as sp

def verify_separatrix():
    """
    This function verifies that the curve d = -u**2 is a separatrix (an invariant manifold)
    for the given system of differential equations.
    """
    print("Verifying the proposed separatrix equation: d = -u**2\n")

    # Define t as a symbol and u as a function of t
    t = sp.Symbol('t')
    u = sp.Function('u')(t)

    # 1. Define the system of ODEs
    # We are given u'(t)
    u_prime_expr = u**2 * (u - 1)
    # And d'(t), where d is also a function of t
    d = sp.Function('d')(t)
    d_prime_expr = 2*d**2 + (5*u**2 - 3*u)*d - u**3*(1-u)

    # 2. Define the proposed separatrix equation
    separatrix_d_expr = -u**2

    # 3. We check if a trajectory following d(t) = -u(t)**2 is consistent with the system.
    #    We check if the d'(t) equation holds true under this assumption.

    # LHS of the equation for d' is the time derivative of d(t)
    # Given d(t) = -u(t)**2, we calculate its derivative with respect to t.
    lhs = sp.diff(separatrix_d_expr, t)
    # The result is -2*u(t)*u'(t). Now substitute the expression for u'(t) from the system.
    lhs_subs = lhs.subs(sp.Derivative(u, t), u_prime_expr)
    
    # RHS of the equation for d' is the given expression.
    # We substitute d(t) = -u(t)**2 into this expression.
    rhs = d_prime_expr.subs(d, separatrix_d_expr)

    # 4. Simplify both sides and compare
    simplified_lhs = sp.simplify(lhs_subs)
    simplified_rhs = sp.simplify(rhs)

    print("LHS (d'(t) based on the separatrix and u'(t)) simplifies to:")
    print(simplified_lhs)
    print("\nRHS (from the d'(t) equation with d=-u**2) simplifies to:")
    print(simplified_rhs)

    # Check if the simplified expressions are identical
    if sp.simplify(simplified_lhs - simplified_rhs) == 0:
        print("\nVerification successful: LHS equals RHS.")
        print("This confirms that d = -u**2 is an invariant manifold of the system.")
    else:
        print("\nVerification failed: LHS does not equal RHS.")

    # 5. Output the equation with its numerical components
    a = -1
    b = 2
    print("\nThe equation of the separatrix is d = a * u**b.")
    print(f"The values of the parameters are: a = {a}, b = {b}.")
    print(f"So, the final equation for the separatrix is: d = {a}*u**{b}")

# Run the verification and print the results
verify_separatrix()