import sympy

def solve_controller():
    """
    Determines and prints the set of all stabilizing controllers H_2(s)
    for H_1(s) = s / (s^2 - 1).
    """
    # Define the symbolic variable s and the parameter K(s)
    s = sympy.Symbol('s')
    # K is represented as a function of s
    K = sympy.Function('K')(s)

    # Plant transfer function H1(s) = N(s) / D(s)
    N = s
    D = s**2 - 1

    # Solution to the Bezout identity N*X + D*Y = 1
    # X = s, Y = -1
    X = s
    Y = -1

    # Youla-Kucera parametrization for the controller H2(s)
    # H2(s) = (X + D*K) / (Y - N*K)
    H2_numerator = X + D * K
    H2_denominator = Y - N * K

    # Expand and collect terms for pretty printing
    num_poly = sympy.Poly(H2_numerator, s, K)
    den_poly = sympy.Poly(H2_denominator, s, K)

    # The prompt asks to output each number in the final equation.
    # We will print the structure of the transfer function.
    
    # Manually format the numerator to match the requested style
    # (s^2 - 1)*K(s) + s
    num_str = f"({sympy.pretty(D, use_unicode=False)})*K(s) + {sympy.pretty(X, use_unicode=False)}"
    
    # Manually format the denominator
    # -s*K(s) - 1
    den_str = f"{sympy.pretty(-N, use_unicode=False)}*K(s) + {sympy.pretty(Y, use_unicode=False)}"

    print("The set of all stabilizing controllers H_2(s) is given by:")
    print(f"H_2(s) = ({num_str}) / ({den_str})")
    print("\nWhere K(s) is a parameter function.")
    print("\nFor the closed-loop system to be stable, K(s) must be a stable and proper rational function.")
    print("For the controller H_2(s) itself to be proper, K(s) must satisfy additional constraints (e.g., be strictly proper with specific leading coefficients).")

solve_controller()