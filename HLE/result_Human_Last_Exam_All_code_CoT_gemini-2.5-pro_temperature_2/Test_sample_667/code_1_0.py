from fractions import Fraction

def solve_yahtzee_problem():
    """
    Calculates and explains the probability difference between two Yahtzee strategies.
    """
    print("This program calculates the difference in probability between two strategies in a game of Yahtzee.")
    print("The dice are infinitely-sided, with P(rolling i) = 2**(-i).\n")

    # --- Strategy A: Keep the '1', aim for five '1's ---
    print("--- Strategy A: Keep the single '1' ---")
    target_A = 1
    dice_to_convert_A = 4
    p_A = Fraction(1, 2**target_A)
    
    # Probability for a single die to become a '1' in at most two rolls
    prob_per_die_A = p_A + (1 - p_A) * p_A
    
    # Total probability for Strategy A
    total_prob_A = prob_per_die_A ** dice_to_convert_A
    
    print(f"Target number: {target_A}")
    print(f"Number of dice to re-roll: {dice_to_convert_A}")
    print(f"Probability of rolling a '{target_A}' in one roll (p_1): {p_A}")
    print(f"Probability of a re-rolled die becoming a '{target_A}' in at most two rolls: {p_A} + (1-{p_A})*{p_A} = {prob_per_die_A}")
    print(f"Total probability of getting five '{target_A}'s (P_A): ({prob_per_die_A})^{dice_to_convert_A} = {total_prob_A}\n")


    # --- Strategy B: Keep the three '3's, aim for five '3's ---
    print("--- Strategy B: Keep the three '3's ---")
    target_B = 3
    dice_to_convert_B = 2
    p_B = Fraction(1, 2**target_B)
    
    # Probability for a single die to become a '3' in at most two rolls
    prob_per_die_B = p_B + (1 - p_B) * p_B

    # Total probability for Strategy B
    total_prob_B = prob_per_die_B ** dice_to_convert_B

    print(f"Target number: {target_B}")
    print(f"Number of dice to re-roll: {dice_to_convert_B}")
    print(f"Probability of rolling a '{target_B}' in one roll (p_3): {p_B}")
    print(f"Probability of a re-rolled die becoming a '{target_B}' in at most two rolls: {p_B} + (1-{p_B})*{p_B} = {prob_per_die_B}")
    print(f"Total probability of getting five '{target_B}'s (P_B): ({prob_per_die_B})^{dice_to_convert_B} = {total_prob_B}\n")

    # --- Final Difference Calculation ---
    print("--- Difference ---")
    # To compare, we find a common denominator
    common_denominator = total_prob_A.denominator * total_prob_B.denominator
    num_A = total_prob_A.numerator * (common_denominator // total_prob_A.denominator)
    num_B = total_prob_B.numerator * (common_denominator // total_prob_B.denominator)

    # Note: P(A) is larger, so we calculate P(A) - P(B) for a positive difference.
    difference = total_prob_A - total_prob_B
    
    print("The chance of winning by keeping the '1' is higher.")
    print("The final difference is P(A) - P(B).")
    print(f"Final Equation: {total_prob_A} - {total_prob_B} = {difference.numerator}/{difference.denominator}")
    print(f"As a decimal, the difference is approximately: {float(difference):.6f}")

solve_yahtzee_problem()
<<<1071/4096>>>