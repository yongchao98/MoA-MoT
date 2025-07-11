import math
from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).
    """
    # Parameters of the problem
    m = 7  # Number of faces on the die
    n = 10 # Number of desired outcomes

    # The optimal strategy requires at least 2 rolls because m < n.
    # The expected number of rolls E can be expressed as:
    # E = sum_{k=1 to inf} P(at least k rolls needed)
    # P(>=1 roll) = 1, P(>=2 rolls) = 1.
    # So, E = 2 + sum_{k=3 to inf} P(at least k rolls needed)
    
    # Step 1: Calculate the probability that more than 2 rolls are needed.
    # This occurs if the outcome of the first 2 rolls is ambiguous.
    # After 2 rolls, we have m^2 = 49 outcomes. An outcome N (from 0 to 48)
    # is ambiguous if the interval [N/49, (N+1)/49] contains a point d/10.
    # This occurs if N = floor(d * 49 / 10) for d in {1, ..., 9}.
    
    m_squared = m * m
    ambiguous_Ns = set()
    for d in range(1, n):
        N = math.floor(d * m_squared / n)
        ambiguous_Ns.add(N)
    
    num_ambiguous_outcomes = len(ambiguous_Ns)
    
    # Probability of needing more than 2 rolls
    prob_fail_after_2 = Fraction(num_ambiguous_outcomes, m_squared)

    # Step 2: For any ambiguous outcome, calculate the probability of needing another roll.
    # If an outcome N is ambiguous, we roll a 3rd time. The new interval
    # is refined by a factor of m. The ambiguity is resolved unless the boundary point
    # falls into one specific sub-interval out of m.
    # The probability of needing a 4th roll, given a 3rd was needed, is 1/m.
    # This pattern continues for all subsequent rolls.
    prob_fail_continues = Fraction(1, m)

    # Step 3: Sum the probabilities for all rolls from the 3rd onwards.
    # The sum is P(>=3) + P(>=4) + P(>=5) + ...
    # This is a geometric series:
    # prob_fail_after_2 * (1 + prob_fail_continues + prob_fail_continues^2 + ...)
    # The sum of the series is 1 / (1 - prob_fail_continues)
    
    geom_series_sum = 1 / (1 - prob_fail_continues)
    
    # The expected number of rolls beyond the first two is prob_fail_after_2 * geom_series_sum
    expected_additional_rolls = prob_fail_after_2 * geom_series_sum

    # The total expected number of rolls is 2 plus this sum.
    E = Fraction(2) + expected_additional_rolls

    # --- Final Output ---
    print("This problem asks for the minimal expected number of rolls of a 7-sided die to uniformly generate a digit from 0 to 9.")
    print("The optimal strategy's expected number of rolls E can be calculated as follows:\n")
    print(f"Number of die faces (m): {m}")
    print(f"Number of target outcomes (n): {n}\n")
    print("At least 2 rolls are always necessary, since 7 < 10.")
    print("After 2 rolls, there are m^2 = 49 outcomes. We calculate the number of outcomes that are 'ambiguous'.")
    print(f"Number of ambiguous outcomes after 2 rolls: {num_ambiguous_outcomes}")
    print(f"Total outcomes after 2 rolls: {m_squared}")
    print(f"Probability of an ambiguous outcome after 2 rolls: {prob_fail_after_2.numerator}/{prob_fail_after_2.denominator}\n")
    print("If an outcome is ambiguous, the probability of needing yet another roll is 1/m = 1/7.")
    print("The expected number of *additional* rolls after an ambiguity occurs is a geometric series sum: 1 / (1 - 1/7) = 7/6.\n")
    print("The final equation for the total expected rolls E is:")
    print(f"E = 2 + (Probability of ambiguity) * (Expected additional rolls)")
    print(f"E = 2 + ({prob_fail_after_2.numerator}/{prob_fail_after_2.denominator}) * ({geom_series_sum.numerator}/{geom_series_sum.denominator})")
    print(f"E = 2 + {expected_additional_rolls.numerator}/{expected_additional_rolls.denominator}")
    print(f"E = {E.numerator}/{E.denominator}\n")
    print("Final simplified result:")
    print(f"{E.numerator}/{E.denominator}")

solve_dice_problem()