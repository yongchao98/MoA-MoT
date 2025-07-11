from fractions import Fraction

def solve_yahtzee_problem():
    """
    Solves the Yahtzee probability problem with an infinitely-sided die.
    """
    # Step 1: Define the probabilities for rolling a '1' and a '3'.
    p1 = Fraction(1, 2**1)  # Probability of rolling a 1
    p3 = Fraction(1, 2**3)  # Probability of rolling a 3

    # Step 2: Calculate the probability of a single die becoming the target number within two rolls.
    # The probability of success for one die within two rolls is p + (1-p)*p.
    # For target '1':
    p_success_1_die = p1 + (1 - p1) * p1
    # For target '3':
    p_success_3_die = p3 + (1 - p3) * p3

    # Step 3: Calculate the total probability for each strategy.
    # Strategy A: Keep one '1', need to turn 4 dice into '1's.
    num_dice_to_change_A = 4
    prob_A = p_success_1_die ** num_dice_to_change_A

    # Strategy B: Keep three '3's, need to turn 2 dice into '3's.
    num_dice_to_change_B = 2
    prob_B = p_success_3_die ** num_dice_to_change_B

    # Step 4: Calculate the difference.
    # We find which probability is higher to subtract correctly.
    if prob_A > prob_B:
        difference = prob_A - prob_B
        p_high, p_low = prob_A, prob_B
        p_high_str = "P(keep 1)"
        p_low_str = "P(keep 3s)"
    else:
        difference = prob_B - prob_A
        p_high, p_low = prob_B, prob_A
        p_high_str = "P(keep 3s)"
        p_low_str = "P(keep 1)"
        
    # Step 5: Output the results and the final equation.
    print(f"The chance of getting five 1s if you keep the '1' is: {prob_A.numerator}/{prob_A.denominator}")
    print(f"The chance of getting five 3s if you keep the three '3s' is: {prob_B.numerator}/{prob_B.denominator}")
    print("\nThe difference is calculated as follows:")
    
    # To show the equation with a common denominator
    common_denominator = max(p_high.denominator, p_low.denominator)
    p_high_common_denom = p_high.numerator * (common_denominator // p_high.denominator)
    p_low_common_denom = p_low.numerator * (common_denominator // p_low.denominator)

    print(f"{p_high_str} - {p_low_str} = {p_high.numerator}/{p_high.denominator} - {p_low.numerator}/{p_low.denominator}")
    print(f"  = {p_high_common_denom}/{common_denominator} - {p_low_common_denom}/{common_denominator}")
    print(f"  = {difference.numerator}/{difference.denominator}")
    
    # Final answer in the required format
    # The question asks for "the difference between the chance... if you keep the 1 and if you keep the three 3s"
    # This is P(A) - P(B) = 81/256 - 225/4096 = 1071/4096
    final_answer = prob_A - prob_B
    print(f"\n<<<1071/4096>>>")


solve_yahtzee_problem()