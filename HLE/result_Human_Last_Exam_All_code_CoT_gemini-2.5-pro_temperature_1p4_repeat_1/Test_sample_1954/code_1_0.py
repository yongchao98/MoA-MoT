def print_minimax_risk_formula():
    """
    This function prints the formula for the minimax risk for estimating the
    parameter theta of a Binomial distribution based on the problem description.

    The derivation shows that the minimax risk R is:
    R = 1 / (4 * (n + 1)^2)
    """

    # The numbers in the final equation
    numerator = 1
    denominator_coefficient = 4
    summand_in_parentheses = 1
    exponent = 2

    print("The minimax risk for estimating theta is given by the formula:")
    print(f"{numerator} / ({denominator_coefficient} * (n + {summand_in_parentheses})^{exponent})")

if __name__ == '__main__':
    print_minimax_risk_formula()