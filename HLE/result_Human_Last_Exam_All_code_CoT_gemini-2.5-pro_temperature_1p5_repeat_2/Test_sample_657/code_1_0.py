import sympy

def solve_controller():
    """
    This function determines the set of all proper stabilizing controllers H2(s)
    for the system H1(s) = s / (s^2 - 1) using Youla-Kucera parametrization.
    """
    # Define symbolic variables
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Plant transfer function H1(s) = n(s)/d(s)
    n = s
    d = s**2 - 1

    # To perform a coprime factorization over RH_infinity, we divide n(s) and d(s)
    # by a stable polynomial of degree 2. We choose p(s) = (s+a)^2.
    # Let's choose a=2 for our calculation.
    a = 2
    p = (s + a)**2

    # We need to solve the polynomial Bezout identity: n*Xp + d*Yp = p
    # From manual calculation (see explanation), we have a particular solution:
    Xp = 5*s + 4
    Yp = -4

    # The set of all stabilizing controllers H2(s) is parametrized by K(s)
    # H2(s) = (Xp + d*K(s)) / (Yp - n*K(s))
    # This formula is a direct result of the Youla-Kucera parametrization.
    # We substitute Xp, Yp, n, d.
    
    # Numerator of H2(s)
    num_H2 = Xp + d * K
    
    # Denominator of H2(s)
    den_H2 = Yp - n * K
    
    # The final expression for the controller H2(s)
    final_H2 = sympy.ratsimp(num_H2 / den_H2)
    
    # Extract expanded numerator and denominator for clear presentation
    num_final, den_final = sympy.fraction(final_H2.expand())
    num_final_expanded = sympy.collect(num_final, K)
    den_final_expanded = sympy.collect(den_final, K)

    print("The set of all stabilizing controllers is given by the transfer function H_2(s):")
    
    # We print each part of the equation as requested
    print("\n\t\t", num_final_expanded)
    print("H_2(s) = \t----------------------------------------")
    print("\t\t", den_final_expanded)

    print("\nwhere K(s) is the parametrization function.")
    print("For internal stability of the closed-loop system, K(s) must be a stable and proper rational function.")
    print("For the controller H_2(s) itself to be proper, K(s) must satisfy the additional condition: lim (s*K(s)) = -5 as s -> infinity.")


solve_controller()