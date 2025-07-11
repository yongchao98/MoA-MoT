def print_stabilizing_controller():
    """
    This function prints the set of all proper stabilizing controllers H_2(s)
    for the plant H_1(s) = s / (s^2 - 1).
    """

    # Numerator and denominator of the plant H_1(s)
    Np = "s"
    Dp = "s**2 - 1"

    # Polynomials X(s) and Y(s) from the Bezout identity
    X = "14*s + 13"
    Y = "s - 8"

    # The Youla-Kucera parametrization for the controller H_2(s)
    numerator = f"({X}) + K(s)*({Dp})"
    denominator = f"({Y}) - K(s)*({Np})"

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print(f"         {numerator}")
    print(f"H_2(s) = {'-' * len(numerator)}")
    print(f"         {denominator}")
    print("\nwhere K(s) is any stable and strictly proper rational function.")
    print("\nEach number in the final equation:")
    print("X(s) = 14*s + 13")
    print("Y(s) = 1*s - 8")
    print("Np(s) = 1*s")
    print("Dp(s) = 1*s**2 - 1")


print_stabilizing_controller()