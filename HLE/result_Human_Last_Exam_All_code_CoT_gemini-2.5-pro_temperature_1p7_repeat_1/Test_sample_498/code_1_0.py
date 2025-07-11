def find_next_letter():
    """
    Solves the letter sequence puzzle based on the QWERTY keyboard layout.
    """
    # The bottom row of letters on a standard QWERTY keyboard
    keyboard_bottom_row = "ZXCVBNM"
    
    # The given sequence from the user
    sequence = "ZXXCVYBN"

    print("The puzzle sequence is based on the bottom row of a QWERTY keyboard.")
    print(f"The standard keyboard row is: {' -> '.join(list(keyboard_bottom_row))}")
    print(f"The given sequence is: {', '.join(list(sequence))}")
    
    # The underlying pattern ignores the inserted 'X' and 'Y'
    base_pattern_followed = "ZXCVBN"
    print(f"The sequence follows the base pattern: {' -> '.join(list(base_pattern_followed))}")

    # Find the last letter from the base pattern in the keyboard row
    last_letter = base_pattern_followed[-1]
    
    # Find its index to determine the next letter
    index_of_last = keyboard_bottom_row.find(last_letter)
    
    # The next letter is at the subsequent index
    if index_of_last + 1 < len(keyboard_bottom_row):
        next_letter = keyboard_bottom_row[index_of_last + 1]
        print(f"The letter after '{last_letter}' on the keyboard is '{next_letter}'.")
    else:
        print("The pattern has reached the end of the keyboard row.")

find_next_letter()