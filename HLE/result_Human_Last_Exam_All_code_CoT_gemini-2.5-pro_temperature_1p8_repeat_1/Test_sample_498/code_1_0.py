def solve_sequence():
    """
    Solves the letter sequence puzzle based on the QWERTY keyboard layout.
    """
    # The bottom row of a standard QWERTY keyboard.
    qwerty_bottom_row = "ZXCVBNM"
    
    # The sequence provided by the user.
    sequence = "ZXXCVYBN"

    # We assume the main pattern is the keyboard row and the extra letters are noise.
    # Find the last character from the user's sequence that exists in our keyboard row.
    last_char_in_pattern = ''
    for char in sequence:
        if char in qwerty_bottom_row:
            last_char_in_pattern = char
            
    # Find the index of this last character.
    try:
        last_index = qwerty_bottom_row.find(last_char_in_pattern)

        # The next letter is the one at the next index.
        if last_index < len(qwerty_bottom_row) - 1:
            next_letter = qwerty_bottom_row[last_index + 1]
            
            # Print the logic step-by-step as requested.
            print("The underlying pattern is the bottom row of the QWERTY keyboard.")
            print("Keyboard Row Pattern:", ", ".join(list(qwerty_bottom_row)))
            print("Given Sequence:", ", ".join(list(sequence)))
            print(f"The last letter in the sequence that follows the pattern is: '{last_char_in_pattern}'")
            print(f"The letter that comes after '{last_char_in_pattern}' in the keyboard row is '{next_letter}'.")

        else:
            print(f"The pattern ends with '{last_char_in_pattern}', which is the last key in the row.")

    except ValueError:
        print("Could not determine the pattern.")

solve_sequence()