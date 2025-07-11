import sympy

def demonstrate_discontinuity():
    """
    This function demonstrates the discontinuity of the minimal polynomial map.
    It considers a derogatory matrix M0 and a sequence of non-derogatory
    matrices Mk converging to M0. It shows that the minimal polynomial of
    Mk does not converge to the minimal polynomial of M0.
    """
    x = sympy.symbols('x')

    # 1. Define a derogatory matrix M0.
    # M0 has eigenvalues {2, 2, 3}. Because the geometric multiplicity of
    # eigenvalue 2 is 2 (two Jordan blocks), it is derogatory.
    M0 = sympy.Matrix([
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 3]
    ])

    # 2. Compute the minimal polynomial of M0.
    # The minimal polynomial should be (x-2)(x-3) = x^2 - 5x + 6.
    pi_M0 = M0.minpoly(x)

    print("The derogatory matrix M0 is:")
    sympy.pprint(M0)
    print("\nIts minimal polynomial pi_M0 is:")
    sympy.pprint(pi_M0)
    print(f"pi_M0(x) = {pi_M0.as_expr()}")
    print("-" * 30)

    # 3. Define a sequence Mk converging to M0.
    # We add a small term 'epsilon' to make the matrix non-derogatory.
    # As epsilon -> 0, Mk -> M0.
    # Let's use a small value for epsilon, e.g., 1/1000.
    epsilon = sympy.Rational(1, 1000)
    Mk = sympy.Matrix([
        [2, epsilon, 0],
        [0, 2,     0],
        [0, 0,     3]
    ])

    # 4. Compute the minimal polynomial of Mk.
    # For any epsilon > 0, Mk has a single Jordan block for the eigenvalue 2,
    # making it non-derogatory. Its minimal polynomial must be equal to its
    # characteristic polynomial: (x-2)^2 * (x-3) = x^3 - 7x^2 + 16x - 12.
    pi_Mk = Mk.minpoly(x)

    print(f"The nearby non-derogatory matrix Mk (with epsilon = {epsilon}) is:")
    sympy.pprint(Mk)
    print("\nIts minimal polynomial pi_Mk is:")
    sympy.pprint(pi_Mk)
    print(f"pi_Mk(x) = {pi_Mk.as_expr()}")
    print("-" * 30)
    
    print("Conclusion:")
    print("As epsilon approaches 0, Mk converges to M0.")
    print("However, the polynomial pi_Mk does not converge to pi_M0.")
    print("The degrees are different (3 vs 2), and the coefficients do not match.")
    print("This demonstrates that the map M -> pi_M is discontinuous at the derogatory matrix M0.")

demonstrate_discontinuity()