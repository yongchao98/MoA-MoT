def solve_sequence():
    """
    Solves the letter sequence puzzle by identifying the underlying pattern.
    """
    # The primary pattern is the bottom row of a QWERTY keyboard.
    qwerty_pattern = "ZXCVBNM"
    
    # The given sequence from the user.
    user_sequence = "ZXXCVYBN"
    
    # The last letter of the given sequence is 'N'.
    last_letter = user_sequence[-1]
    
    try:
        # Find the position of 'N' in the keyboard pattern.
        last_letter_index = qwerty_pattern.index(last_letter)
        
        # The next letter is the one that follows 'N' in the pattern.
        if last_letter_index + 1 < len(qwerty_pattern):
            next_letter = qwerty_pattern[last_letter_index + 1]
            
            # Extract the core pattern from the user's sequence for the final output.
            core_pattern = []
            for letter in qwerty_pattern:
                if letter in user_sequence:
                    core_pattern.append(letter)
            
            print("The underlying pattern is the QWERTY keyboard's bottom row.")
            print("The core sequence is:")
            # Print each letter of the core sequence, emulating an equation as requested.
            equation_str = " -> ".join(core_pattern)
            print(f"{equation_str} -> {next_letter}")
        else:
            print("The pattern is complete.")
            
    except ValueError:
        print(f"Could not determine the next step based on the letter '{last_letter}'.")

solve_sequence()