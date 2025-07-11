import fractions

def solve_yahtzee_problem():
    """
    Solves the Yahtzee probability problem as described.
    Calculates the probabilities for two different strategies and finds their difference.
    """

    print("Analyzing the probability for two strategies in a game of Yahtzee with special dice.\n")

    # --- Strategy A: Keep the '1' ---
    print("--- Strategy A: Keep the '1' and re-roll the other 4 dice ---")
    target_A = 1
    dice_to_roll_A = 4
    
    # Probability of rolling a '1' is p(1) = 1/2^1
    p_A = fractions.Fraction(1, 2**target_A)
    
    # Probability of one die becoming '1' in up to two rolls: 1 - (1 - p)^2
    prob_one_die_A = 1 - (1 - p_A)**2
    
    # Total probability for Strategy A: (prob_one_die_A)^4
    total_prob_A = prob_one_die_A**dice_to_roll_A
    
    print(f"The goal is a Yahtzee of '1's. We need {dice_to_roll_A} dice to become '1'.")
    print(f"The probability of rolling a '1' is p(1) = 1/2.")
    print(f"The probability of one die becoming a '1' in two rolls is 1 - (1 - {p_A})**2 = {prob_one_die_A}.")
    print(f"The total probability for this strategy is ({prob_one_die_A})**{dice_to_roll_A} = {total_prob_A.numerator}/{total_prob_A.denominator}.")
    
    print("\n" + "="*50 + "\n")

    # --- Strategy B: Keep the three '3's ---
    print("--- Strategy B: Keep the three '3's and re-roll the other 2 dice ---")
    target_B = 3
    dice_to_roll_B = 2

    # Probability of rolling a '3' is p(3) = 1/2^3
    p_B = fractions.Fraction(1, 2**target_B)

    # Probability of one die becoming '3' in up to two rolls: 1 - (1 - p)^2
    prob_one_die_B = 1 - (1 - p_B)**2
    
    # Total probability for Strategy B: (prob_one_die_B)^2
    total_prob_B = prob_one_die_B**dice_to_roll_B

    print(f"The goal is a Yahtzee of '3's. We need {dice_to_roll_B} dice to become '3'.")
    print(f"The probability of rolling a '3' is p(3) = {p_B}.")
    print(f"The probability of one die becoming a '3' in two rolls is 1 - (1 - {p_B})**2 = {prob_one_die_B}.")
    print(f"The total probability for this strategy is ({prob_one_die_B})**{dice_to_roll_B} = {total_prob_B.numerator}/{total_prob_B.denominator}.")

    print("\n" + "="*50 + "\n")
    
    # --- Difference Calculation ---
    print("--- Calculating the Difference (P(Strategy A) - P(Strategy B)) ---")
    
    # Find common denominator for display purposes
    common_denominator = max(total_prob_A.denominator, total_prob_B.denominator)
    num_A = total_prob_A.numerator * (common_denominator // total_prob_A.denominator)
    num_B = total_prob_B.numerator * (common_denominator // total_prob_B.denominator)

    difference = total_prob_A - total_prob_B
    
    print(f"Difference = {total_prob_A.numerator}/{total_prob_A.denominator} - {total_prob_B.numerator}/{total_prob_B.denominator}")
    print(f"Difference = {num_A}/{common_denominator} - {num_B}/{common_denominator}")
    print(f"Difference = ({num_A} - {num_B}) / {common_denominator}")
    print(f"Difference = {difference.numerator}/{difference.denominator}")
    
solve_yahtzee_problem()
<<<1071/4096>>>