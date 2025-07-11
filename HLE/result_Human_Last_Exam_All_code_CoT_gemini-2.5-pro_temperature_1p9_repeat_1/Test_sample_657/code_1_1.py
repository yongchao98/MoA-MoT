import sympy

def solve_controller_parametrization():
    """
    This script determines the set of all proper stabilizing controllers H_2(s)
    for the plant H_1(s) = s / (s^2 - 1) using Youla-Kucera parametrization.
    """
    # Define the symbolic variable 's' and the parameter function 'K(s)'
    s = sympy.Symbol('s')
    K_func = sympy.Function('K')
    K = K_func(s)

    # 1. Define the numerator and denominator of the plant H_1(s)
    N = s
    D = s**2 - 1

    # 2. Provide a particular solution to the Bezout identity D(s)X(s) + N(s)Y(s) = 1.
    # For D(s) = s^2 - 1 and N(s) = s, a solution is X(s) = -1, Y(s) = s.
    X = -1
    Y = s

    # 3. Formulate the general expression for the controller H_2(s).
    # H_2(s) = (X(s) + N(s)K(s)) / (Y(s) - D(s)K(s))
    # where K(s) is any stable and proper rational function.
    H2_numerator = X + N * K
    H2_denominator = Y - D * K

    # Collect terms for better readability
    num_expr = sympy.collect(H2_numerator, K)
    den_expr = sympy.collect(H2_denominator, K)

    # 4. Print the final result in the requested format.
    # We print each number and symbol in the final equation.
    # We use sympy's string printer to get the expression.
    num_str = sympy.sstr(num_expr)
    den_str = sympy.sstr(den_expr)

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print(f"H_2(s) = ({num_str}) / ({den_str})")
    print("\nwhere K(s) can be any stable and proper rational function.")


if __name__ == '__main__':
    solve_controller_parametrization()