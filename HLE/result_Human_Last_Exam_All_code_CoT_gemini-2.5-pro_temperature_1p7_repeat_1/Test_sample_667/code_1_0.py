import fractions

def solve_yahtzee_probability():
    """
    Calculates the difference in probabilities between two Yahtzee strategies.
    """

    # Step 1: Define the probability function p(i) = 1/2^i
    # We use the fractions module for exact rational number arithmetic.
    p = lambda i: fractions.Fraction(1, 2**i)

    # Step 2: Analyze Strategy A (Keep the '1')
    # Goal: Get five 1s. We have one, need to convert four dice.
    # Target number is 1.
    target_A = 1
    dice_to_roll_A = 4
    p_A_single_roll = p(target_A)

    # Probability for one die to become a '1' in at most two rolls:
    # P(success) = P(success on 1st) + P(fail on 1st) * P(success on 2nd)
    prob_success_A_one_die = p_A_single_roll + (1 - p_A_single_roll) * p_A_single_roll

    # Total probability for Strategy A (all 4 dice must succeed)
    prob_total_A = prob_success_A_one_die ** dice_to_roll_A

    # Step 3: Analyze Strategy B (Keep the three '3's)
    # Goal: Get five 3s. We have three, need to convert two dice.
    # Target number is 3.
    target_B = 3
    dice_to_roll_B = 2
    p_B_single_roll = p(target_B)

    # Probability for one die to become a '3' in at most two rolls:
    prob_success_B_one_die = p_B_single_roll + (1 - p_B_single_roll) * p_B_single_roll
    
    # Total probability for Strategy B (both 2 dice must succeed)
    prob_total_B = prob_success_B_one_die ** dice_to_roll_B

    # Step 4: Calculate the difference and print the results with the final equation.
    difference = prob_total_A - prob_total_B

    print("--- Strategy A: Keep the '1' ---")
    print(f"Target: five 1s. Probability of rolling a 1 is p(1) = {p_A_single_roll}.")
    print(f"Prob. of one die becoming a '1' in two rolls: {p_A_single_roll} + (1-{p_A_single_roll})*{p_A_single_roll} = {prob_success_A_one_die}")
    print(f"Total probability P(A) = ({prob_success_A_one_die})^4 = {prob_total_A}")
    print("\n--- Strategy B: Keep three '3's ---")
    print(f"Target: five 3s. Probability of rolling a 3 is p(3) = {p_B_single_roll}.")
    print(f"Prob. of one die becoming a '3' in two rolls: {p_B_single_roll} + (1-{p_B_single_roll})*{p_B_single_roll} = {prob_success_B_one_die}")
    print(f"Total probability P(B) = ({prob_success_B_one_die})^2 = {prob_total_B}")
    print("\n--- Difference Calculation ---")
    print(f"The difference is P(A) - P(B).")
    print(f"Difference = {prob_total_A} - {prob_total_B}")
    
    # To show the common denominator calculation
    common_denom = 4096
    prob_A_common_denom = prob_total_A.numerator * (common_denom // prob_total_A.denominator)
    
    print(f"Difference = {prob_A_common_denom}/{common_denom} - {prob_total_B.numerator}/{common_denom}")
    print(f"Difference = {difference.numerator}/{difference.denominator}")
    
    # Final answer in the requested format
    print(f"\n<<<{difference.numerator}/{difference.denominator}>>>")

solve_yahtzee_probability()