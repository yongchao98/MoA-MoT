import math
from fractions import Fraction

def solve_ur_probability():
    """
    Calculates the probability of capturing the opponent's piece in the described scenario.
    """
    # Step 1: Define game parameters based on our plan.
    # The opponent's piece is on our square #5.
    target_square = 5
    # To get an extra turn, we must land on the entry rosette at square #4.
    entry_rosette_square = 4

    # Step 2: Determine the required sequence of dice rolls.
    # First roll must be a 4 to land on the rosette.
    first_roll_needed = entry_rosette_square
    # Second roll must be a 1 to move from square 4 to 5.
    second_roll_needed = target_square - entry_rosette_square

    print("To capture the piece, the following must happen in a single turn:")
    print(f"1. Your first roll must be a {first_roll_needed}. This enters a piece on square {first_roll_needed}, which is a rosette, granting an extra roll.")
    print(f"2. Your second roll must be a {second_roll_needed}, to move your new piece from square {entry_rosette_square} to {target_square} for the capture.\n")

    # Step 3: Calculate the probability of each specific roll.
    # There are 4 tetrahedral dice, so 2^4 = 16 total equally likely outcomes.
    # The probability of rolling a sum 'k' is C(4, k) / 16.
    total_outcomes = 16

    # Probability of the first roll being a 4. (C(4, 4) = 1)
    ways_to_roll_4 = math.comb(4, first_roll_needed)
    prob_first_roll = Fraction(ways_to_roll_4, total_outcomes)

    # Probability of the second roll being a 1. (C(4, 1) = 4)
    ways_to_roll_1 = math.comb(4, second_roll_needed)
    prob_second_roll = Fraction(ways_to_roll_1, total_outcomes)

    print(f"The probability of the first roll being a {first_roll_needed} is {ways_to_roll_4}/{total_outcomes}.")
    print(f"The probability of the second roll being a {second_roll_needed} is {ways_to_roll_1}/{total_outcomes}.\n")

    # Step 4: Calculate the total probability by multiplying the independent event probabilities.
    total_probability = prob_first_roll * prob_second_roll

    print("The final probability is the product of these two probabilities:")
    print(f"{prob_first_roll.numerator}/{prob_first_roll.denominator} * {prob_second_roll.numerator}/{prob_second_roll.denominator} = {total_probability.numerator}/{total_probability.denominator}")
    print(f"\nThe total probability of capturing the piece is: {total_probability}")

solve_ur_probability()