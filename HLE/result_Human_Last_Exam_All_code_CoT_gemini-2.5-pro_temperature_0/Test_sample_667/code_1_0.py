from fractions import Fraction

def solve_yahtzee_probability():
    """
    Calculates and explains the probability difference between two Yahtzee strategies.
    """
    # --- Strategy A: Keep the '1' ---
    # We need to turn 4 dice into '1's. We have two re-rolls.
    target_face_A = 1
    num_dice_to_roll_A = 4
    
    # Probability of rolling a '1' on a single die.
    p_A_target = Fraction(1, 2**target_face_A)
    
    # Probability of a single die becoming a '1' within two rolls.
    p_success_one_die_A = 1 - (1 - p_A_target)**2
    
    # Probability of all 4 dice becoming '1's.
    prob_A = p_success_one_die_A ** num_dice_to_roll_A

    # --- Strategy B: Keep the three '3's ---
    # We need to turn 2 dice into '3's. We have two re-rolls.
    target_face_B = 3
    num_dice_to_roll_B = 2

    # Probability of rolling a '3' on a single die.
    p_B_target = Fraction(1, 2**target_face_B)

    # Probability of a single die becoming a '3' within two rolls.
    p_success_one_die_B = 1 - (1 - p_B_target)**2

    # Probability of both 2 dice becoming '3's.
    prob_B = p_success_one_die_B ** num_dice_to_roll_B

    # --- Difference ---
    difference = prob_A - prob_B

    # --- Output ---
    print("This problem asks for the difference in probabilities of achieving a Yahtzee between two strategies.")
    
    print("\n--- Strategy A: Keep the single '1' ---")
    print("We keep the '1' and re-roll the other 4 dice, aiming for five '1's.")
    print(f"The probability of rolling a '1' is p(1) = 1/{2**target_face_A}.")
    print(f"The probability of one die becoming a '1' within two rolls is 1 - (1 - {p_A_target})**2 = {p_success_one_die_A}.")
    print(f"The total probability of success for this strategy is ({p_success_one_die_A})^{num_dice_to_roll_A} = {prob_A}.")

    print("\n--- Strategy B: Keep the three '3's ---")
    print("We keep the three '3's and re-roll the other 2 dice, aiming for five '3's.")
    print(f"The probability of rolling a '3' is p(3) = 1/{2**target_face_B}.")
    print(f"The probability of one die becoming a '3' within two rolls is 1 - (1 - {p_B_target})**2 = {p_success_one_die_B}.")
    print(f"The total probability of success for this strategy is ({p_success_one_die_B})^{num_dice_to_roll_B} = {prob_B}.")

    print("\n--- Difference in Probabilities ---")
    print(f"The difference is the probability of Strategy A minus the probability of Strategy B.")
    
    # To show the common denominator calculation
    common_denominator = difference.denominator
    prob_A_common_denom = Fraction(prob_A.numerator * (common_denominator // prob_A.denominator), common_denominator)
    
    print(f"Difference = {prob_A} - {prob_B}")
    print(f"Difference = {prob_A_common_denom} - {prob_B}")
    print(f"Difference = {difference}")
    
    print("\nFinal Equation:")
    print(f"({p_success_one_die_A})^{num_dice_to_roll_A} - ({p_success_one_die_B})^{num_dice_to_roll_B} = {prob_A} - {prob_B} = {difference}")

if __name__ == '__main__':
    solve_yahtzee_probability()
<<<1071/4096>>>