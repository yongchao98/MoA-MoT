import math

def solve_ur_probability():
    """
    Calculates the probability of capturing an opponent's piece at the end of the
    middle path in the Royal Game of Ur, given the player has no pieces on the board.
    """

    # --- Game State and Rules Definition ---

    # The player's path consists of 4 private squares followed by the 8 shared squares.
    # The opponent's piece is on the 8th square of the shared middle path.
    # So, the target square is at position 4 (private) + 8 (shared) = 12.
    target_square = 12

    # The game uses four binary dice, each can result in a 0 or a 1.
    num_dice = 4
    max_possible_roll = num_dice * 1

    # --- Analysis of the Required Move ---

    # To capture the piece, the player must land on the target square.
    # Since the player has no pieces on the board, they must enter a new piece.
    # A new piece enters at the square number equal to the dice roll.
    # Therefore, the required dice roll must be equal to the target square's position.
    required_roll = target_square

    print("Analyzing the Royal Game of Ur scenario:")
    print(f"1. The opponent's piece is on the 8th square of the middle path, which is square {target_square} on your total path.")
    print("2. You have no pieces on the board, so your only move is to enter a new piece.")
    print(f"3. To land on square {target_square}, you would need to roll a {required_roll}.")
    print(f"4. The maximum possible roll with {num_dice} binary dice is {max_possible_roll}.")
    print(f"5. Since the required roll ({required_roll}) is greater than the maximum possible roll ({max_possible_roll}), a capture is impossible on this turn.")

    # --- Probability Calculation ---

    # The number of outcomes where the player rolls the required number.
    successful_outcomes = 0

    # The total number of possible outcomes is 2 for each of the 4 dice.
    total_outcomes = 2**num_dice

    # The probability is the ratio of successful outcomes to total outcomes.
    # The result is simplified to a common fraction.
    if successful_outcomes == 0:
        numerator = 0
        denominator = 1
    else:
        # This part of the code is not reached in this scenario but is included for completeness.
        common_divisor = math.gcd(successful_outcomes, total_outcomes)
        numerator = successful_outcomes // common_divisor
        denominator = total_outcomes // common_divisor

    print("\nCalculating the final probability:")
    print(f"Number of ways to roll a {required_roll}: {successful_outcomes}")
    print(f"Total number of possible dice outcomes: {total_outcomes}")
    print("\nThe probability equation is:")
    print(f"P(capture) = (Successful Outcomes) / (Total Outcomes) = {successful_outcomes} / {total_outcomes}")
    print(f"\nThe simplified common fraction is: {numerator}/{denominator}")

solve_ur_probability()