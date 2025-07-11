def solve_sequence():
    """
    Solves the letter sequence puzzle by comparing it to the bottom row of a QWERTY keyboard.
    """
    # The sequence of letters provided by the user.
    sequence = "ZXXCVYBN"

    # The reference sequence is the bottom row of a standard QWERTY keyboard.
    keyboard_bottom_row = "ZXCVBNM"

    print("Analyzing the sequence:", sequence)
    print("The underlying pattern appears to be the bottom row of a QWERTY keyboard.")
    print("Keyboard reference:", keyboard_bottom_row)
    print("-" * 50)
    print("Verifying the pattern step-by-step:")

    # Create an iterator for the keyboard row to track the expected letter.
    keyboard_iter = iter(keyboard_bottom_row)
    
    # This will hold the next character we expect to see from the keyboard pattern.
    next_expected_char_in_pattern = next(keyboard_iter)

    # We iterate through the provided sequence to see how it matches the keyboard row.
    for char_in_sequence in sequence:
        # Check if the character from the input sequence matches the expected pattern character.
        if char_in_sequence == next_expected_char_in_pattern:
            print(f"- '{char_in_sequence}' : Matches the expected pattern character '{next_expected_char_in_pattern}'.")
            # If it matches, we look for the next character in the keyboard row.
            try:
                next_expected_char_in_pattern = next(keyboard_iter)
            except StopIteration:
                # This case is for when the sequence finishes the entire keyboard row.
                next_expected_char_in_pattern = "[End of Reference]"
        else:
            # If it doesn't match, it's an extra character not part of the primary pattern.
            print(f"- '{char_in_sequence}' : Does not match expected '{next_expected_char_in_pattern}'. This is an extra letter.")
            # The expected character remains the same for the next iteration.

    print("-" * 50)
    print("After processing the entire sequence, the next character in the keyboard pattern is:")
    
    # The final answer is the character we were expecting next from the pattern.
    final_answer = next_expected_char_in_pattern
    print(final_answer)

solve_sequence()