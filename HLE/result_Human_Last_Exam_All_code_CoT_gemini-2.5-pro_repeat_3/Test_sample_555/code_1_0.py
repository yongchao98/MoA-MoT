import math
from fractions import Fraction

def solve():
    """
    This function calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).
    """
    min_expected_value = None
    optimal_n = 0

    # A single roll (n=1) has 7 outcomes, which is less than 10, so it's not possible
    # to generate 10 uniform results. We must use n >= 2 rolls per trial.
    # We can analyze E(n) for n from 2 up to a reasonable limit (e.g., 10).
    # E(n) is generally greater than n, so we expect the minimum to occur for a small n.
    for n in range(2, 11):
        # Total number of outcomes for n rolls
        total_outcomes = 7**n
        
        # We need to map outcomes to 10 digits uniformly.
        # Find the largest multiple of 10 less than or equal to total_outcomes.
        usable_outcomes = 10 * math.floor(total_outcomes / 10)
        
        if usable_outcomes == 0:
            continue
            
        # The probability of a successful trial (one that yields a digit)
        # p_success = usable_outcomes / total_outcomes
        
        # The expected number of trials is 1/p_success. Each trial takes n rolls.
        # E(n) = n * (1 / p_success) = n * total_outcomes / usable_outcomes
        expected_value = Fraction(n * total_outcomes, usable_outcomes)
        
        # Keep track of the minimum expected value and the n that produces it.
        if min_expected_value is None or expected_value < min_expected_value:
            min_expected_value = expected_value
            optimal_n = n

    # After finding the optimal n, we print the detailed calculation for that n.
    n = optimal_n
    total_outcomes = 7**n
    floor_value = math.floor(total_outcomes / 10)
    usable_outcomes = 10 * floor_value
    
    unsimplified_numerator = n * total_outcomes
    unsimplified_denominator = usable_outcomes
    
    final_fraction = min_expected_value
    
    print(f"To find the minimal expected value, we analyze the strategy of taking n rolls at a time.")
    print(f"The expected number of rolls is given by the formula E(n) = (n * 7^n) / (10 * floor(7^n/10)).")
    print(f"By testing values of n, we find the minimum is at n = {optimal_n}.")
    print(f"\nFor n = {optimal_n}, the calculation is:")
    print(f"E({optimal_n}) = ({optimal_n} * 7^{optimal_n}) / (10 * floor(7^{optimal_n}/10))")
    print(f"E({optimal_n}) = ({optimal_n} * {total_outcomes}) / (10 * {floor_value})")
    print(f"E({optimal_n}) = {unsimplified_numerator} / {unsimplified_denominator}")
    print(f"The minimal expected value as a simplified fraction is: {final_fraction.numerator}/{final_fraction.denominator}")

solve()
<<<49/20>>>