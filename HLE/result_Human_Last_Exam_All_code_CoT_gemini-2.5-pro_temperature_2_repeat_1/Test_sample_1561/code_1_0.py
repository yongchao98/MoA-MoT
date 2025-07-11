import math

def calculate_capture_probability():
    """
    Calculates the probability of capturing an opponent's piece
    at the end of the middle path on the first turn.
    """
    print("This script calculates the probability of capturing the opponent's piece.")
    
    # Step 1: Explain the interpretation of the game state.
    print("\n--- Step 1: Interpreting the Board ---")
    print("The opponent is at the 'very end of the middle path'.")
    print("In the standard rules, the final square of the middle path (square 12) is a safe rosette.")
    print("Assuming this is a solvable puzzle, we interpret the opponent's position as the last *capturable* square in the middle path, which is square 11.")

    # Step 2: Explain the required sequence of rolls.
    print("\n--- Step 2: Determining the Required Rolls ---")
    print("To capture on square 11 from the start, a specific sequence of rolls leveraging rosettes is needed:")
    print("1. Roll a 4: Land on square 4 (rosette) -> Get an extra roll.")
    print("2. Roll a 4: Move from square 4 to 8 (rosette) -> Get an extra roll.")
    print("3. Roll a 3: Move from square 8 to 11 (capture).")
    
    # Step 3: Calculate the probabilities of individual rolls.
    print("\n--- Step 3: Calculating Individual Roll Probabilities ---")
    num_dice = 4
    total_outcomes = 2**num_dice
    
    # Calculate combinations for rolling a 3 and a 4
    ways_to_roll_3 = math.comb(num_dice, 3)
    ways_to_roll_4 = math.comb(num_dice, 4)
    
    print(f"The game uses {num_dice} binary dice, giving {total_outcomes} possible outcomes per roll.")
    print(f"Probability of rolling a 4 (C(4,4) ways): {ways_to_roll_4}/{total_outcomes}")
    print(f"Probability of rolling a 3 (C(4,3) ways): {ways_to_roll_3}/{total_outcomes}")
    
    # Step 4: Calculate the combined probability.
    print("\n--- Step 4: Calculating the Final Probability ---")
    print("The final probability is the product of the sequence P(roll=4) * P(roll=4) * P(roll=3).")
    
    # Numerator of the final fraction
    final_numerator = ways_to_roll_4 * ways_to_roll_4 * ways_to_roll_3
    # Denominator of the final fraction
    final_denominator = total_outcomes * total_outcomes * total_outcomes

    print(f"Calculation: ({ways_to_roll_4}/{total_outcomes}) * ({ways_to_roll_4}/{total_outcomes}) * ({ways_to_roll_3}/{total_outcomes}) = {final_numerator}/{final_denominator}")

    # Step 5: Simplify the fraction.
    common_divisor = math.gcd(final_numerator, final_denominator)
    simplified_num = final_numerator // common_divisor
    simplified_den = final_denominator // common_divisor
    
    print("\n--- Final Answer ---")
    print(f"The simplified probability of capturing the piece is:")
    print(f"{simplified_num}/{simplified_den}")

calculate_capture_probability()