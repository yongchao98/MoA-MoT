def solve_sequence():
    """
    Solves the letter sequence puzzle by identifying the pattern
    based on a QWERTY keyboard layout.
    """
    # The bottom row of a standard QWERTY keyboard
    qwerty_bottom_row = "ZXCVBNM"

    # The given sequence
    sequence = "ZXXCVYBN"

    # The core pattern identified from the sequence
    core_pattern = "ZXCVBN"

    print(f"The given sequence is: {sequence}_")
    print(f"The primary pattern is based on the bottom row of a QWERTY keyboard: {qwerty_bottom_row}")
    print(f"The sequence follows the keyboard pattern '{core_pattern}', with 'X' and 'Y' as extra letters.")
    
    # Find the last letter of the core pattern on the keyboard row
    last_letter = core_pattern[-1]
    try:
        last_letter_index = qwerty_bottom_row.index(last_letter)
        
        # The next letter is the one immediately following it on the keyboard row
        if last_letter_index + 1 < len(qwerty_bottom_row):
            next_letter = qwerty_bottom_row[last_letter_index + 1]
            print(f"Following this pattern, the letter after '{last_letter}' is '{next_letter}'.")
        else:
            print("The pattern has reached the end of the keyboard row.")

    except ValueError:
        print("The core pattern is not on the keyboard row.")

solve_sequence()