from fractions import Fraction

def find_P_formula():
    """
    This function calculates the coefficients for the correction factor P(n)
    and prints the resulting formula.
    """

    # The Euler-Maclaurin expansion of ln(Q(n)) gives the series for ln(Q(n)/T(n)).
    # The term with n^-2 comes from the B_4 Bernoulli number term.
    # d_2 = B_4 / (4!) * f'''(n) / (-n^-2) = (-1/30)/24 * (-1) = 1/720
    d2_coeff = Fraction(1, 720)

    # The term with n^-4 comes from the B_6 Bernoulli number term.
    # d_4 = B_6 / (6!) * f^(5)(n) / (-n^-4) = (1/42)/720 * (-6) = -1/5040
    d4_coeff = Fraction(-1, 5040)

    # We need P(n) to be exp(d2/n^2 + d4/n^4 + ...).
    # We expand the exponential up to the n^-4 term.
    # P(n) approx 1 + (d2/n^2 + d4/n^4) + (1/2)*(d2/n^2)^2
    # The coefficient of n^-2 in P(n) is d2.
    a2_coeff = d2_coeff

    # The coefficient of n^-4 in P(n) is d4 + d2^2 / 2.
    a4_coeff = d4_coeff + d2_coeff**2 / 2

    # Extract numerators and denominators for clear printing.
    a2_num = a2_coeff.numerator
    a2_den = a2_coeff.denominator

    a4_num = a4_coeff.numerator
    a4_den = a4_coeff.denominator

    # Print the formula with each number clearly visible.
    print("The formula for P(n) is:")
    print(f"P(n) = 1 + {a2_num} / ({a2_den} * n^2) - {abs(a4_num)} / ({a4_den} * n^4)")

if __name__ == '__main__':
    find_P_formula()
