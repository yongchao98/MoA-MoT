def find_next_letter():
    """
    This function solves the sequence puzzle by assuming it follows the
    bottom row of a standard QWERTY keyboard.
    """
    # The full bottom row of a QWERTY keyboard
    keyboard_row = "ZXCVBNM"
    
    # The core pattern recognized from the user's sequence
    sequence_base = "ZXCVBN"
    
    print(f"The sequence is based on the keyboard row: {', '.join(list(keyboard_row))}")
    print(f"The core sequence identified is: {', '.join(list(sequence_base))}")
    
    # Find the last letter from our identified sequence
    last_letter = sequence_base[-1]
    
    # Find the index of this letter in the full keyboard row
    try:
        index_of_last = keyboard_row.find(last_letter)
        
        # Determine the next letter if it exists
        if index_of_last < len(keyboard_row) - 1:
            next_letter = keyboard_row[index_of_last + 1]
            print(f"The letter after '{last_letter}' in the row is '{next_letter}'.")
        else:
            print("The sequence is already at the end of the keyboard row.")
            
    except ValueError:
        print("A letter in the sequence was not found on the keyboard row.")

find_next_letter()