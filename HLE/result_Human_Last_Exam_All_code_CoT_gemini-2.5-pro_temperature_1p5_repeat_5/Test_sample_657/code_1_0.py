import sympy

def solve():
    """
    Determines the set of all proper stabilizing controllers H2(s) for the system H1(s) = s/(s^2 - 1).
    """
    s = sympy.Symbol('s')

    # According to the Youla-Kucera parametrization, the set of all stabilizing controllers
    # H2(s) = U(s)/V(s) is given by the general solution to the Diophantine equation:
    # (s^2 - 1)*V(s) + s*U(s) = A(s)
    # where A(s) is an arbitrary Hurwitz polynomial (all roots have negative real parts).
    #
    # The general solution for U(s) and V(s) can be parametrized by an arbitrary
    # polynomial K(s). A particular solution to s*U_0 + (s^2-1)*V_0 = 1 is U_0=s, V_0=-1.
    # So the general solution for U(s) and V(s) is:
    # U(s) = s*A(s) - K(s)*(s^2 - 1)
    # V(s) = -A(s) + K(s)*s
    #
    # Thus, the controller is H2(s) = U(s)/V(s).

    # Define symbolic functions for K(s) and A(s) for display purposes.
    K = sympy.Function('K')(s)
    A = sympy.Function('A')(s)

    # Numerator of the controller
    U = s * A - K * (s**2 - 1)
    
    # Denominator of the controller
    V = -A + s * K
    
    # Final expression for the controller H2(s)
    H2 = U / V
    
    # --- Final Output ---
    # The problem asks for the answer in the form of a transfer function.
    # We will print the numerator and denominator separately to form the equation.
    
    print("The set of all stabilizing controllers for H1(s) is given by H2(s) = U(s)/V(s), where:")
    
    # We use sympy.pretty to format the polynomials nicely.
    U_pretty = sympy.pretty(U, use_unicode=True)
    V_pretty = sympy.pretty(V, use_unicode=True)

    # Find the length of the horizontal bar for the fraction
    max_len = max(len(U_pretty), len(V_pretty))
    bar = "\u2014" * max_len

    print("\n          " + U_pretty)
    print("H2(s) = " + bar)
    print("          " + V_pretty)

    print("\nHere, A(s) is an arbitrary Hurwitz polynomial (all its roots have negative real parts),")
    print("and K(s) is an arbitrary polynomial.")
    print("\nFor H2(s) to be a *proper* controller, the degrees of the polynomials A(s) and K(s)")
    print("must satisfy constraints that are generally not fulfilled if K(s) is an arbitrary polynomial.")
    print("Properness is only achieved for specific choices of A(s) and K(s) where their degrees and")
    print("coefficients are tightly related, effectively removing K(s) as a free polynomial parameter.")

solve()
<<<H_2(s) = (s*A(s) - K(s)*(s**2 - 1)) / (-A(s) + s*K(s))>>>