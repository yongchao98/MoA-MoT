def solve_sequence():
    """
    Solves the letter sequence puzzle by identifying it as the bottom
    row of a QWERTY keyboard.
    """
    # The letters on the bottom row of a QWERTY keyboard form the pattern.
    keyboard_bottom_row = "ZXCVBNM"
    
    # The input sequence from the problem.
    sequence = "ZXXCVYBN"
    
    # We can filter the input sequence to find only the characters 
    # that match our keyboard row pattern. This helps us find our place.
    # A list comprehension creates a new list containing only valid keys.
    main_pattern_in_sequence = [char for char in sequence if char in keyboard_bottom_row]
    
    # The last letter of this filtered sequence is our reference point.
    last_letter_in_pattern = main_pattern_in_sequence[-1]
    
    # Find the index of this last letter in the complete keyboard row string.
    last_letter_index = keyboard_bottom_row.find(last_letter_in_pattern)
    
    # The next letter in the puzzle is the one at the subsequent index.
    next_letter = keyboard_bottom_row[last_letter_index + 1]
    
    print("The underlying pattern follows the QWERTY keyboard's bottom row.")
    print("The complete pattern can be shown as the following sequence of letters:")
    
    # Output each letter in the final pattern "equation"
    for i, letter in enumerate(keyboard_bottom_row):
        print(letter, end="")
        if i < len(keyboard_bottom_row) - 1:
            print(" -> ", end="")
    print("\n")
    
    print(f"The last letter from the input that fits this pattern is '{last_letter_in_pattern}'.")
    print(f"Therefore, the next letter in the sequence should be '{next_letter}'.")

solve_sequence()