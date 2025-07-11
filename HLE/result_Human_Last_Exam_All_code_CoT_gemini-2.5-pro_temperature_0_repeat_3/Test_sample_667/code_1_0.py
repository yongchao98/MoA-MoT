from fractions import Fraction

def solve_yahtzee_problem():
    """
    Calculates the difference in probabilities for two Yahtzee strategies.
    """
    # --- Strategy A: Keep the single '1' ---
    # We need to roll four '1's. We have two re-rolls.
    # Probability of rolling a '1' on a single die.
    p1 = Fraction(1, 2**1)
    
    # Probability of a single die becoming a '1' within two rolls.
    # This is 1 - (probability of not getting a '1' on either roll).
    prob_one_die_becomes_1 = 1 - (1 - p1)**2
    
    # Probability of all four dice becoming '1's.
    prob_A = prob_one_die_becomes_1**4

    # --- Strategy B: Keep the three '3's ---
    # We need to roll two '3's. We have two re-rolls.
    # Probability of rolling a '3' on a single die.
    p3 = Fraction(1, 2**3)

    # Probability of a single die becoming a '3' within two rolls.
    prob_one_die_becomes_3 = 1 - (1 - p3)**2

    # Probability of both two dice becoming '3's.
    prob_B = prob_one_die_becomes_3**2

    # --- Calculate the difference ---
    difference = prob_A - prob_B

    # --- Print the results ---
    print("This script calculates the difference between two Yahtzee strategies.")
    print("\nStrategy A: Keep the single '1'.")
    print(f"The probability of one die becoming a '1' in two rolls is: {prob_one_die_becomes_1.numerator}/{prob_one_die_becomes_1.denominator}")
    print(f"The total probability of success P(A) is ({prob_one_die_becomes_1.numerator}/{prob_one_die_becomes_1.denominator})^4 = {prob_A.numerator}/{prob_A.denominator}")

    print("\nStrategy B: Keep the three '3's.")
    print(f"The probability of one die becoming a '3' in two rolls is: {prob_one_die_becomes_3.numerator}/{prob_one_die_becomes_3.denominator}")
    print(f"The total probability of success P(B) is ({prob_one_die_becomes_3.numerator}/{prob_one_die_becomes_3.denominator})^2 = {prob_B.numerator}/{prob_B.denominator}")

    print("\n--- Final Calculation ---")
    print(f"The difference is P(A) - P(B):")
    print(f"{prob_A.numerator}/{prob_A.denominator} - {prob_B.numerator}/{prob_B.denominator} = {difference.numerator}/{difference.denominator}")
    print(f"\nAs a decimal, the difference is: {float(difference)}")

solve_yahtzee_problem()
<<<1071/4096>>>