import math
from fractions import Fraction

def solve():
    """
    This function calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).
    """
    # The strategy is rejection sampling in blocks. We analyze the expected
    # number of rolls E(n) for a block size of n.
    # E(n) = n * (7^n) / (10 * floor(7^n / 10))
    # We must have n >= 2, since 7^1 < 10.

    # It can be proven that E(n) is an increasing function for n >= 2.
    # Therefore, the minimum expected value occurs at n=2.
    n = 2

    rolls_per_try = n
    total_outcomes = 7**n
    used_outcomes = 10 * (total_outcomes // 10)

    # The expected number of rolls is E(2) = 2 * (49 / 40)
    unsimplified_numerator = rolls_per_try * total_outcomes
    unsimplified_denominator = used_outcomes

    # Using the Fraction class to find the simplified form
    simplified_fraction = Fraction(unsimplified_numerator, unsimplified_denominator)
    final_numerator = simplified_fraction.numerator
    final_denominator = simplified_fraction.denominator

    print("The optimal strategy uses rejection sampling with a block size of n=2 rolls.")
    print("Here is the calculation for the expected number of rolls:")
    print(f"E(2) = {rolls_per_try} * (7^{n} / (10 * floor(7^{n}/10)))")
    print(f"     = {rolls_per_try} * ({total_outcomes} / {used_outcomes})")
    print(f"     = {unsimplified_numerator} / {unsimplified_denominator}")
    print(f"     = {final_numerator} / {final_denominator}")
    print("\nThe minimal expected value of rolls is expressed as the simplified fraction above.")

solve()
<<<49/20>>>