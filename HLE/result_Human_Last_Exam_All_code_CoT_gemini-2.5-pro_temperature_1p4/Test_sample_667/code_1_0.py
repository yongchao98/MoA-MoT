from fractions import Fraction

def solve_yahtzee_probability():
    """
    Calculates and prints the probability difference between two Yahtzee strategies.
    """

    # Step 1: Define the probability function for the die
    def p(i):
        return Fraction(1, 2**i)

    # Step 2: Define a function for the probability of a single die becoming a target
    # in two rolls. P(success) = p_i * (2 - p_i)
    def prob_success_two_rolls(target_i):
        pi = p(target_i)
        return pi * (2 - pi)

    # Step 3: Analyze Strategy A: Keep the single '1' and go for a Yahtzee of 1s.
    target_A = 1
    num_dice_to_roll_A = 4
    prob_one_die_A = prob_success_two_rolls(target_A)
    prob_A = prob_one_die_A ** num_dice_to_roll_A
    p_A = p(target_A)

    print("--- Strategy A: Keep the single '1' ---")
    print(f"The goal is a Yahtzee of 1s. We must turn 4 dice into '1's.")
    print(f"The probability of rolling a 1 is p(1) = 1/{2**target_A} = {p_A}.")
    print(f"The probability for one die to become a '1' in two rolls is:")
    print(f"  p(1) * (2 - p(1)) = {p_A} * (2 - {p_A}) = {prob_one_die_A}")
    print(f"The total probability for Strategy A is ({prob_one_die_A})^4 = {prob_A}.")
    print()

    # Step 4: Analyze Strategy B: Keep the three '3's and go for a Yahtzee of 3s.
    target_B = 3
    num_dice_to_roll_B = 2
    prob_one_die_B = prob_success_two_rolls(target_B)
    prob_B = prob_one_die_B ** num_dice_to_roll_B
    p_B_val = p(target_B)

    print("--- Strategy B: Keep the three '3's ---")
    print(f"The goal is a Yahtzee of 3s. We must turn 2 dice into '3's.")
    print(f"The probability of rolling a 3 is p(3) = 1/{2**target_B} = {p_B_val}.")
    print(f"The probability for one die to become a '3' in two rolls is:")
    print(f"  p(3) * (2 - p(3)) = {p_B_val} * (2 - {p_B_val}) = {prob_one_die_B}")
    print(f"The total probability for Strategy B is ({prob_one_die_B})^2 = {prob_B}.")
    print()

    # Step 5: Calculate and display the difference.
    # Note: Keeping the '1' gives a higher chance, so we calculate P(A) - P(B).
    difference = prob_A - prob_B
    
    # To show the common denominator calculation
    common_denominator = difference.denominator
    prob_A_common_denom = prob_A.limit_denominator(common_denominator)
    prob_B_common_denom = prob_B.limit_denominator(common_denominator)

    print("--- Difference in Probabilities ---")
    print(f"The difference is P(Strategy A) - P(Strategy B):")
    print(f"  {prob_A} - {prob_B}")
    print(f"= {prob_A_common_denom} - {prob_B_common_denom}")
    print(f"= {difference}")
    print(f"\nThe numeric value of the difference is approximately {float(difference):.4f}.")

solve_yahtzee_probability()
<<<1071/4096>>>