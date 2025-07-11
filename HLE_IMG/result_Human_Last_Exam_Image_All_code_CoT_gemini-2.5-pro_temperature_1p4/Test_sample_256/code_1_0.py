def solve_rock_art_puzzle():
    """
    Analyzes the rock art image to identify any symbols not belonging to
    ancient Southwest cultures.
    """

    # Based on visual inspection, the statement is true.
    is_true = True

    # The identified symbol is modern text graffiti.
    symbol = "the word 'no'"

    # The location of the symbol in the image.
    location = "in the upper right quadrant"

    # Construct the final answer.
    if is_true:
        # We format a clear and descriptive answer.
        answer = f"True. {location.capitalize()} of the image, it is possible to distinguish {symbol}. This appears to be modern text and is not a religious symbol from the ancient Southwest cultures that created the primary artwork."
    else:
        # This part of the logic is not reached for this specific problem.
        answer = "False."

    print(answer)

solve_rock_art_puzzle()