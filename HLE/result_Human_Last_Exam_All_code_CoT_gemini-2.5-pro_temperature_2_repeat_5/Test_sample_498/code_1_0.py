def find_next_letter():
    """
    This function identifies the next letter in the sequence based on the
    QWERTY keyboard layout.
    """
    # The bottom row of letters on a standard QWERTY keyboard
    keyboard_bottom_row = "ZXCVBNM"
    
    # The recognizable part of the input sequence ends with 'N'
    last_known_letter = 'N'
    
    print(f"The sequence is based on the bottom letter row of a QWERTY keyboard.")
    print(f"The full sequence of letters is: Z, X, C, V, B, N, M")
    print(f"The given sequence's pattern ends with the letter: {last_known_letter}")
    
    try:
        # Find the index of the last known letter
        index = keyboard_bottom_row.find(last_known_letter)
        
        # Ensure the letter is in the sequence and not the last one
        if index != -1 and index < len(keyboard_bottom_row) - 1:
            # The next letter is at the next index
            next_letter = keyboard_bottom_row[index + 1]
            print(f"The letter that comes after '{last_known_letter}' in this sequence is '{next_letter}'.")
        else:
            print("The end of the known sequence was reached or the letter wasn't found.")

    except Exception as e:
        print(f"An error occurred: {e}")

find_next_letter()
