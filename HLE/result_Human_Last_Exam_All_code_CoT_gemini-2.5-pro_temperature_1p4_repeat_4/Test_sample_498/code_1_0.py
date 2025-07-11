def solve_keyboard_sequence():
    """
    Identifies the next letter in the sequence ZXXCVYBN_.
    The logic is based on the bottom row of a QWERTY keyboard.
    """
    
    # The bottom letter row on a standard QWERTY keyboard.
    qwerty_bottom_row = "ZXCVBNM"
    
    # The core letters of the sequence follow this pattern.
    # The original sequence is "ZXXCVYBN". The extra 'X' and the 'Y' are ignored.
    core_pattern = "Z -> X -> C -> V -> B -> N"
    
    # The last letter in the established pattern is 'N'.
    last_letter = 'N'
    
    # Find the index of the last letter in the keyboard row.
    last_letter_index = qwerty_bottom_row.find(last_letter)
    
    # The next letter is the one at the next position.
    next_letter = qwerty_bottom_row[last_letter_index + 1]

    print("The primary pattern follows the letters on the bottom row of a QWERTY keyboard.")
    print(f"The pattern is: {core_pattern} -> ?")
    print(f"The next letter in this sequence is '{next_letter}'.")

solve_keyboard_sequence()