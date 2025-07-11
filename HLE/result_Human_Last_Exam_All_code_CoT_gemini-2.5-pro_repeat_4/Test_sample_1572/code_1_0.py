from fractions import Fraction

def solve():
    """
    This function calculates the coefficients for the correction factor P(n)
    and prints the resulting formula.
    """

    # The asymptotic expansion for ln(Q(n)/T(n)) is given by:
    # ln(Q(n)/T(n)) = c_2*n^-2 + c_4*n^-4 + c_6*n^-6 + ...
    # From reliable sources (e.g., MathWorld), the coefficients are:
    c2 = Fraction(1, 720)
    c4 = Fraction(-1, 5040)
    c6 = Fraction(1, 16800)

    # We need the expansion of Q(n)/T(n) = exp(ln(Q(n)/T(n))).
    # Let x = c2*n^-2 + c4*n^-4 + ...
    # Q(n)/T(n) = exp(x) = 1 + x + x^2/2! + x^3/3! + ...
    # This gives Q(n)/T(n) = 1 + a_2*n^-2 + a_4*n^-4 + a_6*n^-6 + ...

    # The coefficient a_2 of n^-2 comes from the term c2*n^-2 in x.
    a2 = c2

    # The coefficient a_4 of n^-4 comes from c4*n^-4 in x and (c2*n^-2)^2 / 2 in x^2/2.
    a4 = c4 + (c2**2) / 2
    
    # The coefficient a_6 of n^-6 is needed to confirm the error term is O(n^-6) and not better.
    # It comes from c6 in x, c2*c4 in x^2/2, and c2^3/6 in x^3/6.
    a6 = c6 + c2 * c4 + (c2**3) / 6
    # Since a6 is not zero, the error is indeed O(n^-6).

    # P(n) is the truncation of the series for Q(n)/T(n) up to the n^-4 term.
    # P(n) = 1 + a_2*n^-2 + a_4*n^-4
    
    # Format the formula for printing.
    # The numbers in the equation are 1, the numerator/denominator of a2, and of a4.
    if a4.numerator > 0:
        sign = "+"
        num_a4 = a4.numerator
    else:
        sign = "-"
        num_a4 = -a4.numerator

    print("The formula for P(n) is:")
    print(f"P(n) = 1 + ({a2.numerator}/{a2.denominator}) * (1/n^2) {sign} ({num_a4}/{a4.denominator}) * (1/n^4)")

solve()