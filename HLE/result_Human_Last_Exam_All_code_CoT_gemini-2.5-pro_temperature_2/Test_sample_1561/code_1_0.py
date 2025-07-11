import math
from fractions import Fraction

def get_roll_probability(roll_value, num_dice=4):
    """
    Calculates the probability of a specific roll in the Royal Game of Ur.
    The game uses 4 tetrahedral dice, each with 2 marked and 2 unmarked sides.
    The roll value is the number of marked sides facing up.
    
    Args:
        roll_value (int): The desired roll value (0, 1, 2, 3, or 4).
        num_dice (int): The number of dice thrown (always 4 in this game).
        
    Returns:
        Fraction: The probability of the roll.
    """
    # Total possible outcomes = 2 states (marked/unmarked) ^ num_dice
    total_outcomes = 2**num_dice
    
    # Number of successful outcomes is the binomial coefficient "num_dice choose roll_value"
    successful_outcomes = math.comb(num_dice, roll_value)
    
    return Fraction(successful_outcomes, total_outcomes)

def solve_ur_probability():
    """
    Solves the specific probability problem and prints the step-by-step reasoning.
    """
    # Step 1: Define the required rolls for the successful capture sequence.
    # To capture the piece on square 5 in one turn from off the board, we must:
    # - Roll a 4 to land on the first rosette (square 4) and get an extra turn.
    # - Roll a 1 on the extra turn to move from square 4 to 5 (4 + 1 = 5).
    first_roll_needed = 4
    second_roll_needed = 1
    
    # Step 2: Calculate the probability of each required roll.
    prob_first_roll = get_roll_probability(first_roll_needed)
    prob_second_roll = get_roll_probability(second_roll_needed)
    
    # Step 3: Calculate the total probability by multiplying the independent event probabilities.
    total_probability = prob_first_roll * prob_second_roll

    # Step 4: Print the detailed explanation and the result.
    print("To capture the opponent's piece at square 5 on your turn, a specific sequence is required:")
    print(f"1. You must first roll a '{first_roll_needed}' to enter a piece on the rosette at square 4, earning an extra turn.")
    print(f"2. On your extra turn, you must roll a '{second_roll_needed}' to move the piece from square 4 to 5.")
    print("\nLet's calculate the probabilities:")
    
    # Print the probability calculation for the first roll
    p1_num = prob_first_roll.numerator
    p1_den = prob_first_roll.denominator
    print(f"The probability of rolling a {first_roll_needed} is {p1_num}/{p1_den}.")

    # Print the probability calculation for the second roll
    p2_num = prob_second_roll.numerator
    p2_den = prob_second_roll.denominator
    print(f"The probability of rolling a {second_roll_needed} is {p2_num}/{p2_den}.")

    print("\nThe total probability is the product of these two probabilities.")
    print(f"Final Probability = P(Roll is {first_roll_needed}) * P(Roll is {second_roll_needed})")
    print(f"Final Probability = {p1_num}/{p1_den} * {p2_num}/{p2_den} = {total_probability.numerator}/{total_probability.denominator}")

# Run the solver
solve_ur_probability()
