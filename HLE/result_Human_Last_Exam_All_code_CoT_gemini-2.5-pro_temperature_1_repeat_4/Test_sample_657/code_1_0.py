import sympy

def find_stabilizing_controllers():
    """
    Determines the set of all proper stabilizing controllers H_2(s) for the plant
    H_1(s) = s / (s^2 - 1) using the Youla-Kucera parametrization.
    """
    # Define the symbolic variable 's' and the Youla parameter K(s)
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # The plant H_1(s) = s / (s^2 - 1) has an unstable pole at s = 1.

    # Step 1: Find a coprime factorization H_1(s) = N(s) / D(s) over the ring
    # of stable, proper rational functions (RH_infinity).
    # We select a factorization that isolates the unstable part.
    # N(s) = s / (s + 1)^2
    # D(s) = (s - 1) / (s + 1)
    # Both N(s) and D(s) are stable and proper.
    # N(s)/D(s) = [s/(s+1)^2] / [(s-1)/(s+1)] = s/((s+1)(s-1)) = H_1(s).

    # Step 2: Find a particular solution (X(s), Y(s)) in RH_infinity for the
    # Bezout identity: N(s)*X(s) + D(s)*Y(s) = 1.
    # By solving the Diophantine equation, a simple particular solution is found:
    X_s = 4
    Y_s = (s - 1) / (s + 1)
    # Both X(s) and Y(s) are stable and proper.

    # Step 3: The set of all stabilizing controllers H_2(s) is parametrized by K(s),
    # where K(s) is any stable and proper rational function.
    # H_2(s) = (X(s) + D(s)*K(s)) / (Y(s) - N(s)*K(s))
    
    # We substitute the expressions for N, D, X, Y.
    # Numerator of H2 is 4 + [(s-1)/(s+1)]*K(s)
    # Denominator of H2 is [(s-1)/(s+1)] - [s/(s+1)^2]*K(s)

    # To simplify, we can multiply the numerator and denominator by (s+1)^2.
    num_H2 = 4 * (s + 1)**2 + (s - 1) * (s + 1) * K
    den_H2 = (s - 1) * (s + 1) - s * K

    # Expand and simplify the expressions for clarity.
    final_num = sympy.expand(num_H2)
    final_den = sympy.expand(den_H2)

    print("The set of all proper stabilizing controllers H_2(s) is given by a transfer function parametrized by K(s),")
    print("where K(s) is any stable and proper rational function.\n")
    print("H_2(s) = Numerator / Denominator\n")
    
    print("Numerator:")
    sympy.pprint(final_num)
    print("\nDenominator:")
    sympy.pprint(final_den)

    # To explicitly output each number in the final equation, we can express H_2(s) as:
    # H_2(s) = (P1(s) + P2(s)*K(s)) / (P3(s) + P4(s)*K(s))
    # and print the coefficients of the polynomials P1, P2, P3, P4.
    
    P1 = sympy.poly(4 * (s + 1)**2, s)
    P2 = sympy.poly((s - 1) * (s + 1), s)
    P3 = sympy.poly((s - 1) * (s + 1), s)
    P4 = sympy.poly(-s, s)

    print("\nTo show all numbers in the equation, here are the coefficients of the constituent polynomials:")
    print(f"Numerator term P1(s) = {P1.expr}, Coefficients: {P1.all_coeffs()}")
    print(f"Numerator term P2(s) = {P2.expr}, Coefficients: {P2.all_coeffs()}")
    print(f"Denominator term P3(s) = {P3.expr}, Coefficients: {P3.all_coeffs()}")
    print(f"Denominator term P4(s) = {P4.expr}, Coefficients: {P4.all_coeffs()}")


find_stabilizing_controllers()