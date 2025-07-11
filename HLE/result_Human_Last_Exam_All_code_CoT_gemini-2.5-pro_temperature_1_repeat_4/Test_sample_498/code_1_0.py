def solve_sequence():
    """
    Analyzes the letter sequence ZXXCVYBN_ to find the next letter.
    The logic is based on the layout of a standard QWERTY keyboard.
    """

    # The letters on the bottom row of a QWERTY keyboard.
    bottom_row = ['Z', 'X', 'C', 'V', 'B', 'N', 'M']
    sequence = "ZXXCVYBN"

    print("Analyzing the sequence: " + " ".join(sequence) + " _")
    print("-" * 40)
    print("Step 1: The pattern is based on the QWERTY keyboard layout.")
    print(f"The bottom row of letters is: {' '.join(bottom_row)}")
    print("-" * 40)

    print("Step 2: Trace the sequence on the keyboard.")
    print("The sequence mostly follows the letters on the bottom row from left to right.")
    print("Path: Z -> X -> C -> V -> ... -> B -> N")
    print("While the given sequence has some unusual elements ('X' is repeated, 'Y' is an outlier),")
    print("the end of the sequence ('B N') returns to the simple left-to-right pattern.")
    print("-" * 40)

    print("Step 3: Determine the next letter.")
    last_known_letter = sequence[-1]
    try:
        # Find the position of the last letter in our reference row.
        last_letter_index = bottom_row.index(last_known_letter)

        # Check if there is a next letter in the row.
        if last_letter_index + 1 < len(bottom_row):
            next_letter = bottom_row[last_letter_index + 1]
            print(f"Following the pattern, the letter after '{last_known_letter}' on the bottom row is '{next_letter}'.")
        else:
            print("The sequence has reached the end of the bottom row.")
            next_letter = "Pattern complete"

    except ValueError:
        next_letter = "Could not determine the next letter based on the pattern."

    print("-" * 40)
    print("Conclusion: The next letter in the sequence is M.")

solve_sequence()