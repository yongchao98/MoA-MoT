import math

def calculate_capture_probability():
    """
    Calculates the probability of capturing an opponent's piece in the Royal Game of Ur
    under the specified conditions.
    """
    # Game parameters based on Irving L. Finkel's rules
    num_dice = 4
    # Each die is binary, with outcomes 0 or 1
    max_roll_per_die = 1
    
    # The opponent's piece is at the end of the 8-square shared middle path,
    # which is the 12th square of a player's journey.
    target_square = 12
    
    # To capture the piece, we must land on it. Since we have no pieces on the board,
    # we must enter a new piece. The entry square is determined by the dice roll.
    # Therefore, the required roll is equal to the target square number.
    required_roll = target_square
    
    # The maximum possible roll is the sum of the max value of all dice.
    max_possible_roll = num_dice * max_roll_per_die
    
    # The total number of possible outcomes for the dice roll is 2^4.
    total_outcomes = 2**num_dice
    
    # Determine the number of ways to achieve the required roll.
    # Since the required roll is impossible, the number of successful outcomes is 0.
    if required_roll > max_possible_roll:
        successful_outcomes = 0
    else:
        # This case is not met here, but included for completeness.
        # The number of ways to get sum 'k' with 'n' binary dice is "n choose k".
        successful_outcomes = math.comb(num_dice, required_roll)

    # The probability is the ratio of successful outcomes to total outcomes.
    numerator = successful_outcomes
    denominator = total_outcomes
    
    # Simplify the fraction for the final answer.
    if numerator == 0:
        simplified_denominator = 1
    else:
        common_divisor = math.gcd(numerator, denominator)
        simplified_denominator = denominator // common_divisor

    print("To capture the piece on square 12, you must roll a 12.")
    print(f"The maximum possible roll with {num_dice} binary dice is {max_possible_roll}.")
    print("Since a roll of 12 is impossible, the number of successful outcomes is 0.")
    print(f"The total number of possible outcomes is {total_outcomes}.")
    print("\nThe probability is calculated as: Successful Outcomes / Total Outcomes")
    print(f"Final Equation: P(capture) = {numerator} / {denominator}")
    print(f"Simplified Fraction: {numerator}/{simplified_denominator}")

calculate_capture_probability()
<<<0/1>>>