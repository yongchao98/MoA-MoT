import fractions

def solve_yahtzee_problem():
    """
    Calculates the difference in probabilities for two strategies in a game of Yahtzee
    with infinitely-sided dice.
    """

    # --- Scenario A: Keep the three 3s ---

    # Probability of rolling a 3
    p3 = fractions.Fraction(1, 2**3)

    # We need to roll two dice and have them both become 3s. We have two attempts.
    # For a single die, the probability of it becoming a 3 in two rolls is:
    # P(success on 1st) + P(fail on 1st) * P(success on 2nd)
    prob_one_die_becomes_3 = p3 + (1 - p3) * p3

    # Since we need this to happen for two independent dice, we square the probability.
    prob_A = prob_one_die_becomes_3 ** 2

    # --- Scenario B: Keep the single 1 ---

    # Probability of rolling a 1
    p1 = fractions.Fraction(1, 2**1)

    # We need to roll four dice and have them all become 1s. We have two attempts.
    # For a single die, the probability of it becoming a 1 in two rolls is:
    prob_one_die_becomes_1 = p1 + (1 - p1) * p1

    # Since we need this to happen for four independent dice, we raise the probability to the 4th power.
    prob_B = prob_one_die_becomes_1 ** 4

    # --- Calculate the difference ---
    # The question asks for the difference between the two chances.
    difference = prob_B - prob_A

    # --- Print the results ---
    print("This script calculates the difference between two Yahtzee strategies.")
    print("\nScenario A: Keep the three 3s.")
    print(f"The probability of success is {prob_A.numerator}/{prob_A.denominator}")

    print("\nScenario B: Keep the single 1.")
    print(f"The probability of success is {prob_B.numerator}/{prob_B.denominator}")

    # To show the final equation, we find a common denominator for the two probabilities.
    common_denominator = prob_A.denominator
    # prob_B needs to be converted to the same denominator for the printout.
    prob_B_common_denom_numerator = prob_B.numerator * (common_denominator // prob_B.denominator)

    print("\nThe difference between the chance of getting Yahtzee if you keep the 1 and if you keep the three 3s is:")
    print(f"{prob_B} - {prob_A} = {prob_B_common_denom_numerator}/{common_denominator} - {prob_A.numerator}/{common_denominator} = {difference.numerator}/{difference.denominator}")
    print(f"The final difference is {difference.numerator}/{difference.denominator}")

solve_yahtzee_problem()
<<<1071/4096>>>