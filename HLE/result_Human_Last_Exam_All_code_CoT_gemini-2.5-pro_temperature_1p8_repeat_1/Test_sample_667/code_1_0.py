from fractions import Fraction

def calculate_probability():
    """
    Calculates the difference in probabilities between two Yahtzee strategies.
    """

    # p(i) is the probability of rolling 'i' in a single throw
    # p(i) = 2**(-i)

    # For any die, the probability of it becoming the target number 'i' in
    # at most two rolls is P_success(i) = p(i) + (1-p(i))*p(i).
    def prob_one_die_success(i):
        p_i = Fraction(1, 2**i)
        return p_i + (1 - p_i) * p_i

    # Strategy 1: Keep the '1'. We need the other four dice to become '1's.
    # The probability for a single die to become a '1' is:
    prob_s1_per_die = prob_one_die_success(1)
    # The total probability for this strategy is this value to the power of 4.
    prob_s1 = prob_s1_per_die ** 4

    # Strategy 2: Keep the three '3's. We need the other two dice to become '3's.
    # The probability for a single die to become a '3' is:
    prob_s2_per_die = prob_one_die_success(3)
    # The total probability for this strategy is this value to the power of 2.
    prob_s2 = prob_s2_per_die ** 2

    # Calculate the difference
    difference = prob_s1 - prob_s2

    print("--- Yahtzee Probability Analysis ---")
    print("\nStrategy 1: Keep the single '1' and roll for four more '1's.")
    print(f"The probability of one die becoming '1' in two rolls is {prob_s1_per_die.numerator}/{prob_s1_per_die.denominator}")
    print(f"The total chance of success is ({prob_s1_per_die.numerator}/{prob_s1_per_die.denominator})^4 = {prob_s1.numerator}/{prob_s1.denominator}")

    print("\nStrategy 2: Keep the three '3's and roll for two more '3's.")
    print(f"The probability of one die becoming '3' in two rolls is {prob_s2_per_die.numerator}/{prob_s2_per_die.denominator}")
    print(f"The total chance of success is ({prob_s2_per_die.numerator}/{prob_s2_per_die.denominator})^2 = {prob_s2.numerator}/{prob_s2.denominator}")

    print("\n--- Final Calculation ---")
    print(f"The difference between the probabilities is:")
    print(f"P(Strategy 1) - P(Strategy 2) = {prob_s1.numerator}/{prob_s1.denominator} - {prob_s2.numerator}/{prob_s2.denominator}")
    
    # To show the common denominator calculation
    common_denominator = max(prob_s1.denominator, prob_s2.denominator)
    s1_num_common = prob_s1.numerator * (common_denominator // prob_s1.denominator)
    s2_num_common = prob_s2.numerator * (common_denominator // prob_s2.denominator)
    
    print(f"                           = {s1_num_common}/{common_denominator} - {s2_num_common}/{common_denominator}")
    print(f"                           = {difference.numerator}/{difference.denominator}")
    print(f"                           ≈ {float(difference):.5f}")

calculate_probability()