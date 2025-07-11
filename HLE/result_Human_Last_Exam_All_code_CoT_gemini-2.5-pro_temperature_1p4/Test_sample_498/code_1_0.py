def solve_sequence():
    """
    Solves the letter sequence puzzle by identifying the underlying pattern
    from the QWERTY keyboard layout.
    """
    # The sequence is based on the bottom row of a QWERTY keyboard.
    bottom_row = "ZXCVBNM"
    
    # The given sequence with some distracting letters.
    input_sequence = "ZXXCVYBN"
    
    # We will extract the letters that follow the correct order of the bottom row.
    # This filters out the extra 'X' and the 'Y'.
    extracted_pattern = ""
    last_index = -1
    for char in input_sequence:
        try:
            current_index = bottom_row.index(char)
            # Only add the character if it appears after the previous one in the sequence.
            if current_index > last_index:
                extracted_pattern += char
                last_index = current_index
        except ValueError:
            # Character is not in the bottom row (like 'Y'), so we ignore it.
            continue
            
    # The last character of our filtered pattern is the last valid step.
    last_char_in_pattern = extracted_pattern[-1]
    
    # Find the index of this last character in the full bottom row.
    final_index = bottom_row.index(last_char_in_pattern)
    
    # The next letter is at the next index.
    if final_index + 1 < len(bottom_row):
        next_letter = bottom_row[final_index + 1]
        print(f"The extracted pattern is: {extracted_pattern}")
        print(f"The last letter in the pattern is '{last_char_in_pattern}'.")
        print(f"The next letter on the keyboard's bottom row is: {next_letter}")
    else:
        print("The sequence is already at the end of the keyboard row.")

solve_sequence()