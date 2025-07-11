def solve_caesar_riddle():
    """
    Calculates the length of the longest possible message Caesar can write.
    """
    # The maximum number of characters the encrypted message can have.
    max_encrypted_length = 10000

    # To get the longest plaintext message, we need the most efficient encoding.
    # The most efficient encoding uses the shortest possible code for each character.
    # The absolute shortest code to represent one character is another single character.
    # This means the minimum possible length for an encoded character is 1.
    # For example: 'E' -> 'I'
    min_code_length_per_char = 1

    # The length of the longest message is the maximum encrypted length
    # divided by the minimum code length per character.
    longest_message_length = max_encrypted_length / min_code_length_per_char

    # Output the explanation and the final equation.
    print(f"The paper has a capacity for {max_encrypted_length} encrypted characters.")
    print(f"To maximize the original message's length, each character must be encoded using the shortest possible sequence of Roman numerals.")
    print(f"The minimum possible code length for a single character is {min_code_length_per_char}.")
    print("\nThe final equation is: Longest Message Length = Max Encrypted Length / Min Code Length")
    # Print the equation with the numbers plugged in.
    print(f"Longest Message Length = {max_encrypted_length} / {min_code_length_per_char}")

    print(f"\nTherefore, the length of the longest message is {int(longest_message_length)}.")

solve_caesar_riddle()
<<<10000>>>