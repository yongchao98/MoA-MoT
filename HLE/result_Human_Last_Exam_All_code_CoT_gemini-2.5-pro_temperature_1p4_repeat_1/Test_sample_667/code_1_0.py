from fractions import Fraction

def calculate_yahtzee_probability():
    """
    Calculates the difference in probabilities between two Yahtzee strategies.
    """

    # --- Scenario A: Keep the single '1' ---
    # We need to turn 4 dice into '1's.
    # The probability of rolling a '1' is p(1) = 2**(-1)
    p1 = Fraction(1, 2)
    # The probability of a single die becoming a '1' in two rolls is q1 = p1 + (1-p1)*p1
    q1 = p1 + (1 - p1) * p1
    # We need this to happen for 4 independent dice.
    prob_A = q1 ** 4

    # --- Scenario B: Keep the three '3's ---
    # We need to turn 2 dice into '3's.
    # The probability of rolling a '3' is p(3) = 2**(-3)
    p3 = Fraction(1, 8)
    # The probability of a single die becoming a '3' in two rolls is q3 = p3 + (1-p3)*p3
    q3 = p3 + (1 - p3) * p3
    # We need this to happen for 2 independent dice.
    prob_B = q3 ** 2

    # --- Calculate the difference ---
    # To subtract, find a common denominator. 4096 is the common denominator.
    prob_A_common_denom = Fraction(prob_A.numerator * (prob_B.denominator // prob_A.denominator), prob_B.denominator)
    difference = prob_A - prob_B
    
    print("This script calculates the difference in probabilities for two strategies in a game of Yahtzee with infinitely-sided dice.")
    print("\n--- Strategy 1: Keep the '1' ---")
    print(f"The probability of one die becoming a '1' in two rolls is: {q1.numerator}/{q1.denominator}")
    print(f"The total probability of success is ({q1.numerator}/{q1.denominator})^4 = {prob_A.numerator}/{prob_A.denominator}")

    print("\n--- Strategy 2: Keep the three '3's ---")
    print(f"The probability of one die becoming a '3' in two rolls is: {q3.numerator}/{q3.denominator}")
    print(f"The total probability of success is ({q3.numerator}/{q3.denominator})^2 = {prob_B.numerator}/{prob_B.denominator}")
    
    print("\n--- The Difference ---")
    print("The equation for the difference is:")
    print(f"{prob_A.numerator}/{prob_A.denominator} - {prob_B.numerator}/{prob_B.denominator} = {prob_A_common_denom.numerator}/{prob_A_common_denom.denominator} - {prob_B.numerator}/{prob_B.denominator} = {difference.numerator}/{difference.denominator}")
    print(f"\nThe numerical value of the difference is: {float(difference)}")

calculate_yahtzee_probability()
<<<1071/4096>>>