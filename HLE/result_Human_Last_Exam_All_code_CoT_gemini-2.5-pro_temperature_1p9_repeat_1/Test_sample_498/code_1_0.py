def solve_keyboard_sequence():
    """
    Determines the next letter in the sequence based on the QWERTY keyboard layout.
    """
    # The bottom letter row of a standard QWERTY keyboard.
    keyboard_row = "ZXCVBNM"
    
    # The sequence ZXXCVYBN_ largely follows the keyboard row.
    # The last key from this pattern in the sequence is 'N'.
    last_key_in_sequence = 'N'
    
    # Find the position of 'N' in the keyboard row.
    position_of_last_key = keyboard_row.find(last_key_in_sequence)
    
    # The next letter is at the next position in the row.
    next_letter = keyboard_row[position_of_last_key + 1]
    
    print("The primary pattern is the bottom row of a QWERTY keyboard:")
    print("Z -> X -> C -> V -> B -> N -> ?")
    print(f"The last key in the sequence from this pattern is '{last_key_in_sequence}'.")
    print(f"The key after '{last_key_in_sequence}' on the keyboard is '{next_letter}'.")
    
    # Displaying the logic as a final equation
    step_value = 1
    print("\nFinal Equation:")
    print(f"{last_key_in_sequence} + {step_value} key position = {next_letter}")

solve_keyboard_sequence()