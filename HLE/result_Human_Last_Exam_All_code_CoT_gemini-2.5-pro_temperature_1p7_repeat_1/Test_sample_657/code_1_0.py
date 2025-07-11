def solve_controller():
    """
    This function determines and prints the set of all proper stabilizing controllers
    for the plant H1(s) = s / (s^2 - 1).
    The method used is the Youla-Kucera parametrization.
    """

    # Numerator of the plant H1(s)
    Np = "s"
    # Denominator of the plant H1(s)
    Dp = "s^2 - 1"

    # A particular solution (X, Y) to the Bezout identity Np*X + Dp*Y = 1
    # s*X + (s^2 - 1)*Y = 1
    # Solution X = s, Y = -1 is found by inspection.
    X = "s"
    Y = "-1"

    # The general form of the stabilizing controller H2(s) is (X + Dp*K) / (Y - Np*K)
    # where K(s) is any stable and strictly proper rational function.

    # Numerator of the controller H2(s)
    num_H2 = f"({X}) + ({Dp})*K(s)"

    # Denominator of the controller H2(s)
    den_H2 = f"({Y}) - ({Np})*K(s)"

    print("The set of all proper stabilizing controllers H_2(s) for the plant H_1(s) = s/(s^2-1) is given by:")
    print("H_2(s) = N(s) / D(s), where K(s) is any stable and strictly proper rational function.")
    
    # We explicitly state the coefficients as requested.
    print("\nThe final formula, with each number specified, is:")
    print(f"N(s) = (1)*s + ((1)*s^2 - (1))*K(s)")
    print(f"D(s) = (-1) - (1)*s*K(s)")


solve_controller()