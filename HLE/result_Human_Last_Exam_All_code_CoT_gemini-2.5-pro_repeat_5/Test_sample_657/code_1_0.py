import sympy

def solve_controller():
    """
    This function determines and prints the set of all proper stabilizing controllers
    for the plant H1(s) = s / (s^2 - 1) using Youla-Kucera parametrization.
    """
    # Define the symbolic variable 's' for the transfer function
    s = sympy.Symbol('s')

    # K(s) is the Youla parameter, which can be any stable and proper rational function.
    # We define it as a symbolic function of s for representation.
    K = sympy.Function('K')(s)

    # Based on the Youla-Kucera parametrization for the given plant,
    # the controller H2(s) = N2(s) / D2(s) is derived.

    # Numerator of the controller H2(s)
    N2 = 4 * (s + 1)**2 + K * (s**2 - 1)

    # Denominator of the controller H2(s)
    D2 = (s**2 - 1) - s * K

    # Expand the polynomials to show all the terms and coefficients clearly.
    N2_expanded = sympy.expand(N2)
    D2_expanded = sympy.expand(D2)

    # --- Output the result ---
    print("The set of all proper stabilizing controllers H2(s) is given by the formula:")
    print("H2(s) = N(s) / D(s)\n")
    print("The numerator is:")
    print("N(s) = ", end="")
    sympy.pprint(N2_expanded, use_unicode=False)
    print("\nThe denominator is:")
    print("D(s) = ", end="")
    sympy.pprint(D2_expanded, use_unicode=False)
    print("\n")
    print("In this expression, K(s) is the parametrization function, which can be any stable and proper rational function.")
    print("The controller H2(s) is guaranteed to be proper for any such choice of K(s).")

if __name__ == '__main__':
    solve_controller()