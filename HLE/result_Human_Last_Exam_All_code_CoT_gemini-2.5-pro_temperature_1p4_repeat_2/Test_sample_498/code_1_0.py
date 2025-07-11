def solve_sequence():
    """
    Solves the letter sequence puzzle based on the QWERTY keyboard layout.
    """
    # The sequence of letters on the bottom row of a standard QWERTY keyboard
    keyboard_bottom_row = "ZXCVBNM"
    
    # The input sequence contains the core pattern with some distractor letters.
    # The core pattern follows the keyboard row: Z, X, C, V, B, N...
    print("The pattern is based on the letters on the bottom row of a standard QWERTY keyboard.")
    print("The provided sequence 'ZXXCVYBN' contains distractor letters, but the underlying pattern is a simple progression.")

    # The last letter of the recognized pattern is 'N'.
    last_letter_in_pattern = 'N'
    
    # Find the index of this last letter in our reference keyboard row.
    try:
        last_letter_index = keyboard_bottom_row.index(last_letter_in_pattern)
        
        # The next letter is the one at the following index.
        next_letter = keyboard_bottom_row[last_letter_index + 1]
        
        # To fulfill the "output each number in the final equation" request, 
        # we will display the logical progression of the sequence.
        
        # Get the part of the sequence we've identified so far.
        identified_sequence = keyboard_bottom_row[:last_letter_index + 1]
        
        print("\nThe logical sequence is:")
        # Print the final equation showing the full progression
        final_equation = " -> ".join(list(identified_sequence)) + f" -> {next_letter}"
        print(final_equation)
        
    except (ValueError, IndexError):
        print("Could not determine the next letter in the sequence.")

solve_sequence()