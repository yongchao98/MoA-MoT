def solve_sequence():
    """
    Determines the next letter in the sequence ZXXCVYBN based on the
    QWERTY keyboard layout.
    """
    # The sequence from the problem.
    input_sequence = "ZXXCVYBN"

    # The pattern source: the bottom row of a standard QWERTY keyboard.
    keyboard_bottom_row = "ZXCVBNM"

    print("Analyzing the sequence:", input_sequence)
    print("The underlying pattern is based on the bottom row of a QWERTY keyboard:", keyboard_bottom_row)
    print("-" * 20)
    
    # We find the last character from the input that is part of the keyboard pattern.
    # In this case, we can see the pattern progresses right up to the end of the input string.
    # The last character of the consistent pattern is 'N'.
    last_char_in_pattern = 'N'
    
    print(f"The main pattern Z-X-C-V-B-N is followed, with some interruptions.")
    print(f"The last letter that follows this keyboard pattern is '{last_char_in_pattern}'.")

    # Find the index of this last character in our definitive keyboard row string.
    try:
        index_of_last_char = keyboard_bottom_row.index(last_char_in_pattern)
    except ValueError:
        print("Error: The last character was not found in the keyboard pattern.")
        return

    # The next letter is the one at the next index, if it exists.
    if index_of_last_char < len(keyboard_bottom_row) - 1:
        next_letter = keyboard_bottom_row[index_of_last_char + 1]
        print(f"The next letter after '{last_char_in_pattern}' on the keyboard row is '{next_letter}'.")
        print("\nTherefore, the next letter in the sequence is:")
        print(next_letter)
    else:
        print("The sequence has reached the end of the keyboard pattern.")

solve_sequence()