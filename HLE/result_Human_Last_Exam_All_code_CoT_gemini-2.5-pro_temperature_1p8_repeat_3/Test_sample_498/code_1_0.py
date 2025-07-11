def solve_sequence():
    """
    Solves the letter sequence puzzle.

    The logic is based on the bottom row of a QWERTY keyboard.
    The sequence 'ZXXCVYBN' is an interleaved sequence.

    1. The primary sequence is the letters of the bottom row of a QWERTY keyboard in order:
       Z, X, C, V, B, N, M
    2. A secondary sequence ('X', 'Y') is interleaved into the primary one.
    3. The puzzle sequence is built like this:
       - Z (Primary)
       - X (Primary)
       - X (Secondary)
       - C (Primary)
       - V (Primary)
       - Y (Secondary)
       - B (Primary)
       - N (Primary)
    4. The next letter should continue the primary sequence. After 'N' comes 'M'.
    """

    # The bottom row of a QWERTY keyboard
    qwerty_bottom_row = "ZXCVBNM"

    # The letters from the primary sequence present in the puzzle
    main_sequence_in_puzzle = ['Z', 'X', 'C', 'V', 'B', 'N']

    # The last letter of the primary sequence found in the puzzle
    last_letter = main_sequence_in_puzzle[-1]

    # Find the index of this last letter in the full keyboard row
    last_letter_index = qwerty_bottom_row.find(last_letter)

    # The next letter is the one at the next index
    if last_letter_index != -1 and last_letter_index < len(qwerty_bottom_row) - 1:
        next_letter = qwerty_bottom_row[last_letter_index + 1]
    else:
        # Should not happen in this case, but good for completeness
        next_letter = "Not found"

    print("The given sequence is: Z, X, X, C, V, Y, B, N, _")
    print(f"The primary pattern is the bottom row of a QWERTY keyboard: {', '.join(list(qwerty_bottom_row))}")
    print(f"The last letter from the primary pattern in the sequence is '{last_letter}'.")
    print(f"The next letter in the QWERTY bottom row sequence is '{next_letter}'.")


solve_sequence()