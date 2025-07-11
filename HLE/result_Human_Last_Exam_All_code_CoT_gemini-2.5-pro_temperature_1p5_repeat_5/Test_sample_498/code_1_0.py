def solve_sequence():
    """
    Solves the letter sequence puzzle based on a QWERTY keyboard layout.
    """
    # The pattern is based on the bottom row of letters on a QWERTY keyboard.
    keyboard_bottom_row = "ZXCVBNM"
    
    # The sequence from the problem.
    sequence = "ZXXCVYBN"
    
    print(f"The pattern is based on the bottom row of a QWERTY keyboard, which is: {keyboard_bottom_row}")
    
    # The last letter in the provided sequence that fits the main pattern is 'N'.
    last_letter_in_sequence = 'N'
    
    # Find the index of this last letter within the keyboard row pattern.
    try:
        last_letter_index = keyboard_bottom_row.index(last_letter_in_sequence)
        
        # The next letter is the one at the subsequent index.
        if last_letter_index + 1 < len(keyboard_bottom_row):
            next_letter = keyboard_bottom_row[last_letter_index + 1]
            
            print(f"The given sequence '{sequence}' follows this keyboard pattern, ending on the letter '{last_letter_in_sequence}'.")
            print(f"The next letter in the keyboard pattern is '{next_letter}'.")
            
            # Displaying the final completed sequence as the "final equation".
            completed_sequence = list(sequence + next_letter)
            print("The final sequence is:")
            print(" -> ".join(completed_sequence))
        else:
            print("The sequence is already at the end of the pattern.")

    except ValueError:
        print(f"The letter '{last_letter_in_sequence}' was not found in the pattern.")

solve_sequence()