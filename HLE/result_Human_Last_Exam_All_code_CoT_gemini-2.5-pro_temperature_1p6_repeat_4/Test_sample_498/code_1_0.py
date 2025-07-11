def solve_sequence():
    """
    Solves the letter sequence puzzle by identifying the pattern
    as the bottom row of a QWERTY keyboard.
    """
    # The bottom letter row on a standard QWERTY keyboard
    bottom_row = "ZXCVBNM"

    # The input sequence from the puzzle
    input_sequence = "ZXXCVYBN"
    
    # The letters from the input that match the keyboard row pattern
    core_pattern_letters = []
    for char in bottom_row:
        if char in input_sequence:
            core_pattern_letters.append(char)

    print(f"The puzzle sequence is: {', '.join(list(input_sequence))}")
    print("This sequence is based on the bottom row of a QWERTY keyboard.")
    print(f"The full bottom row is: {', '.join(list(bottom_row))}")
    
    # The last letter of the recognized pattern is N
    last_letter_in_pattern = 'N'
    print(f"The core pattern letters found are: {', '.join(core_pattern_letters)}")
    print(f"The last letter in this pattern is '{last_letter_in_pattern}'.")

    # Find the index of the last letter in the full keyboard row
    try:
        last_letter_index = bottom_row.find(last_letter_in_pattern)
        
        # Ensure the letter is found and is not the last one on the row
        if last_letter_index < len(bottom_row) - 1:
            next_letter = bottom_row[last_letter_index + 1]
            print(f"The letter that comes after '{last_letter_in_pattern}' on the keyboard is '{next_letter}'.")
            print(f"Therefore, the final sequence is: {', '.join(list(bottom_row))}")
        else:
            # This case handles if 'M' was the last letter
            next_letter = "End of row"
            print("The pattern has reached the end of the keyboard row.")
            
    except ValueError:
        next_letter = "Error: Pattern not recognized"
        print(next_letter)

solve_sequence()