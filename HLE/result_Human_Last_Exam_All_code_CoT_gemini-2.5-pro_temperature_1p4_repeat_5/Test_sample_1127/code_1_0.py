def solve_minimal_polynomial():
    """
    This function prints the minimal polynomial for the connective constant of the specified graph.

    The derivation is based on the known value of the connective constant for the
    diced lattice, mu = sqrt(3 + sqrt(3)).
    From t = mu, we get the equation: t^4 - 6*t^2 + 6 = 0.
    """

    # Coefficients of the minimal polynomial P(t) = a_4*t^4 + a_3*t^3 + a_2*t^2 + a_1*t + a_0
    a4 = 1
    a3 = 0
    a2 = -6
    a1 = 0
    a0 = 6

    print("The minimal polynomial P(t) for the connective constant is defined by the equation P(t) = 0.")
    print("The coefficients of the polynomial are:")
    print(f"  Coefficient of t^4: {a4}")
    print(f"  Coefficient of t^3: {a3}")
    print(f"  Coefficient of t^2: {a2}")
    print(f"  Coefficient of t^1: {a1}")
    print(f"  Constant term:      {a0}")

    print("\nThe minimal polynomial equation is:")
    print("t^4 - 6*t^2 + 6 = 0")

solve_minimal_polynomial()