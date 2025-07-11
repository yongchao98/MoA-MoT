def decode_word(word):
    """
    Decodes a word by shifting each letter one position backward in the alphabet.
    It prints the decoding equation for each letter and the final result.
    """
    decoded_word = []
    final_equation_parts = []
    
    print("Decoding the word:", word)
    print("Pattern: Shift each letter back by 1 (e.g., b -> a, a -> z).")
    print("-" * 20)
    
    for char in word:
        # Convert character to its 1-26 alphabetical position
        original_pos = ord(char) - ord('a') + 1
        
        # Shift back by one
        new_pos = original_pos - 1
        
        # Handle wraparound from 'a' to 'z'
        if new_pos == 0:
            new_pos = 26
            
        # Convert new position back to a character
        new_char = chr(new_pos + ord('a') - 1)
        
        # Store the character for the final word
        decoded_word.append(new_char)
        
        # Print the equation for the current character
        print(f"{char} ({original_pos}) - 1 = {new_char} ({new_pos})")
        
    final_word = "".join(decoded_word)
    print("-" * 20)
    print("Final decoded word:", final_word)
    return final_word

# The word to decode
encoded_word = "anhaeketn"
decode_word(encoded_word)