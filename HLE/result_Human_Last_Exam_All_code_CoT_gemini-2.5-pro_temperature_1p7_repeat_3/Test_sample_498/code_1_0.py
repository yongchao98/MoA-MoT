def solve_keyboard_sequence():
    """
    Analyzes the sequence ZXXCVYBN_ by relating it to a QWERTY keyboard.
    """
    # Define the letter rows of a standard QWERTY keyboard
    top_row =    "Q W E R T Y U I O P"
    middle_row = " A S D F G H J K L "
    bottom_row_full = "Z X C V B N M"

    # The sequence provided by the user
    user_sequence = "ZXXCVYBN"

    print("To solve this, let's look at a standard QWERTY keyboard layout:")
    print(f"Top row:    {top_row}")
    print(f"Middle row: {middle_row}")
    print(f"Bottom row:  {bottom_row_full}")
    print("-" * 40)

    print(f"The user's sequence is: {user_sequence}")
    print("\nIf we ignore the second 'X' and the 'Y', we can see a clear pattern.")
    
    # Extract the base pattern from the user sequence
    base_pattern = "ZXCVBN"
    print(f"The underlying pattern is: {base_pattern}")
    
    print(f"\nThis pattern '{base_pattern}' exactly matches the beginning of the keyboard's bottom row: '{bottom_row_full}'")
    
    # Find the next letter in the keyboard sequence
    try:
        last_letter_index = bottom_row_full.find(base_pattern[-1])
        if last_letter_index != -1 and last_letter_index + 1 < len(bottom_row_full):
            next_letter = bottom_row_full[last_letter_index + 1]
            print(f"\nThe sequence is progressing along the bottom row of the keyboard.")
            print(f"The last identified letter in the pattern is '{base_pattern[-1]}'.")
            print(f"The next letter after '{base_pattern[-1]}' on the keyboard is '{next_letter}'.")
        else:
            # This part should not be reached if the logic is correct
            print("\nCould not determine the next letter from the pattern.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_keyboard_sequence()