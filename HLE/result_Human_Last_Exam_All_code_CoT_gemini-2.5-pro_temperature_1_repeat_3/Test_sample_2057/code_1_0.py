def solve_hat_puzzle():
    """
    This function provides the solution to the hat puzzle.
    The logic is derived step-by-step as explained above.
    """

    # The number of people who say "Yes" in the third round.
    # Based on the reasoning, the 4 people with white hats deduce their hat color.
    num_yes = 4

    # The distribution of hats around the table.
    # This is one possible configuration that fits the puzzle's constraints.
    # Other valid answers are rotations or reflections of this sequence.
    # B = Black, W = White
    hat_distribution = ["B", "W", "B", "W", "B", "W", "B", "B", "W"]

    print(f"How many people replied 'Yes.'?")
    print(f"Answer: {num_yes}")
    print("\nHow are black and white hats distributed around the table?")
    print("Answer: One possible distribution is as follows:")
    print(" ".join(hat_distribution))

solve_hat_puzzle()