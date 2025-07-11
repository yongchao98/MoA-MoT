def decode_and_explain(word):
    """
    Decodes a word by shifting each character one position back in the alphabet ('a' becomes 'z').
    It also prints the equation for each character's transformation.
    """
    decoded_word = ""
    print(f"Decoding the word: '{word}'")
    print("The pattern is a Caesar cipher with a shift of -1.")
    print("The equation for each letter is: new_letter = chr(ord(original_letter) - 1)")
    print("-" * 30)

    final_equation_parts = []

    for char in word:
        original_ord = ord(char)
        if char == 'a':
            # Handle wraparound for 'a'
            decoded_char = 'z'
            new_ord = ord('z')
        else:
            decoded_char = chr(original_ord - 1)
            new_ord = original_ord - 1
        
        # Build the visual equation for each character
        equation_step = f"{char}({original_ord}) - 1 = {decoded_char}({new_ord})"
        print(equation_step)
        
        decoded_word += decoded_char
        final_equation_parts.append(f"{decoded_char}")
        
    print("-" * 30)
    print("Final decoded word:")
    print(decoded_word)
    
    # Per the instruction "output each number in the final equation!", we print the components.
    # We will show the final decoded letters which result from the equations.
    print("\nLetters in the final decoded word from each equation:")
    print(" -> ".join(final_equation_parts))

# The word to be decoded
word_to_decode = "anhaeketn"

# Run the decoding process
decode_and_explain(word_to_decode)