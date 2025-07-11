def solve_sequence():
    """
    Solves the letter sequence puzzle based on the QWERTY keyboard layout.
    """
    # The bottom row of a standard QWERTY keyboard
    bottom_row = "ZXCVBNM"
    
    # The given sequence
    sequence = "ZXXCVYBN"

    # The underlying pattern ignores the extra 'X' and 'Y'
    base_pattern = "Z, X, C, V, B, N"
    
    # The last character of the clear pattern is 'N'
    last_char_in_pattern = 'N'
    
    # Find the index of the last character in the keyboard row
    try:
        last_char_index = bottom_row.index(last_char_in_pattern)
        
        # The next character is the one at the next index
        next_char = bottom_row[last_char_index + 1]

        print(f"The sequence is based on the bottom row of a QWERTY keyboard: {bottom_row}")
        print(f"The primary pattern is the letters in order: {base_pattern}")
        print(f"The last letter from this pattern in the sequence is '{last_char_in_pattern}'.")
        print(f"The next letter on the keyboard after '{last_char_in_pattern}' is '{next_char}'.")

    except (ValueError, IndexError):
        print("Could not determine the next letter based on the keyboard pattern.")

solve_sequence()