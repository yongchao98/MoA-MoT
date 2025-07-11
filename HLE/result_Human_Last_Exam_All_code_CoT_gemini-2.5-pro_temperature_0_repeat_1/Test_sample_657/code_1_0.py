import sympy

def solve_controller():
    """
    This function determines and prints the set of all proper stabilizing controllers H_2(s)
    for the system H_1(s) = s / (s^2 - 1).
    The controller H_2(s) is parametrized by a stable, strictly proper function K(s).
    """
    s = sympy.Symbol('s')
    # K represents the transfer function K(s)
    K = sympy.Function('K')(s)

    # The plant H_1(s) has an unstable pole at s=1.
    # Any stabilizing controller H_2(s) must cancel this pole,
    # leading to a parametrization of the form:
    # H_2(s) = ( (s^2 - 1) * K(s) ) / ( (s + 1) - s * K(s) )
    # where K(s) is any stable and strictly proper rational function.

    # Numerator of H_2(s)
    numerator = (s**2 - 1) * K

    # Denominator of H_2(s)
    denominator = (s + 1) - s * K

    print("The set of all proper stabilizing controllers H_2(s) is given by the transfer function:")
    print(f"H_2(s) = Numerator / Denominator")
    print("\nWhere K(s) is any stable and strictly proper rational function.")
    print("\nNumerator:")
    # The prompt asks to output each number in the final equation.
    # We will print the expression in a structured way.
    num_expr = f"({s**2} - 1)*K(s)"
    print(num_expr)

    print("\nDenominator:")
    den_expr = f"({s} + 1) - {s}*K(s)"
    print(den_expr)

solve_controller()