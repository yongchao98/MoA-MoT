import sympy

def solve_controller_parametrization():
    """
    This function determines and prints the set of all stabilizing controllers
    for the system H1(s) = s / (s^2 - 1) using Youla-Kucera parametrization.
    """
    # Set up the symbolic variable 's' and the parameter function 'K(s)'
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Define the numerator and denominator of the plant H1(s)
    N = s
    D = s**2 - 1

    # Solve the Bezout identity: N*X + D*Y = 1.
    # We use a specific solution pair (X, Y) that facilitates finding a proper controller.
    # This pair is derived from a simpler solution (X=s, Y=-1) by transformation.
    # X_final = s - 1*(s**2 - 1) = -s**2 + s + 1
    # Y_final = -1 + 1*s = s - 1
    X_final = -s**2 + s + 1
    Y_final = s - 1

    # The general formula for the stabilizing controller H_2(s) is:
    # H_2(s) = (X + D*K) / (Y - N*K)
    # where K(s) is any stable and proper rational function.
    num_H2 = X_final + D * K
    den_H2 = Y_final - N * K

    # Rearrange the expressions to clearly show the parametrization by K(s)
    num_expr = sympy.collect(num_H2, K)
    den_expr = sympy.collect(den_H2, K)
    
    # --- Output the results ---
    print("The set of all proper stabilizing controllers H_2(s) is parametrized by a stable, proper rational function K(s).")
    print("The general form of the controller is H_2(s) = Numerator / Denominator, where:")
    print("\n")
    
    # Print the terms with explicit coefficients as requested.
    print("Numerator = (1 * s**2 - 1) * K(s) + (-1 * s**2 + 1 * s + 1)")
    print("\n")
    print("Denominator = (-1 * s) * K(s) + (1 * s - 1)")

    # Print the full symbolic expression for clarity
    print("\nWhich can be written as:")
    final_eq = sympy.Eq(sympy.Symbol('H_2(s)'), num_expr / den_expr)
    sympy.pprint(final_eq, use_unicode=False)


if __name__ == '__main__':
    solve_controller_parametrization()
<<<H_2(s) = (K(s)*(s**2 - 1) - s**2 + s + 1)/(K(s)*(-s) + s - 1)>>>