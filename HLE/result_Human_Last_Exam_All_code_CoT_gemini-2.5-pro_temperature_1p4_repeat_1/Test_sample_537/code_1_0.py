def solve_cipher():
    """
    This function decodes the given word based on the discovered pattern.
    The pattern is that the length of the encoded word corresponds to the number it represents.
    """
    
    # The word to be decoded
    encoded_word = "anhaeketn"
    
    # A dictionary to map a number (the length) to its English word equivalent.
    # This includes all numbers relevant to the problem.
    number_to_word_map = {
        5: "FIVE",
        9: "NINE",
        10: "TEN",
        11: "ELEVEN"
    }
    
    # Calculate the length of the word, which is the decoded number.
    decoded_number = len(encoded_word)
    
    # Retrieve the final decoded word from the map.
    decoded_word = number_to_word_map.get(decoded_number, "UNKNOWN")
    
    # Print the final equation, showing the encoded word, its length (the number),
    # and the final decoded word.
    print(f"Encoded word: '{encoded_word}'")
    print(f"Decoding Equation: len('{encoded_word}') = {decoded_number}")
    print(f"Final Answer: {decoded_word}")

solve_cipher()