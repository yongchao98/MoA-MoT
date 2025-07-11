def solve_crypto():
    """
    Solves the encrypted phrase based on logical deduction.
    The final answer is derived from the following mapping, which was deduced step-by-step.
    """
    # Mapping of hex codes to letters based on the deduction
    # This mapping is the core result of the analysis
    mapping = {
        '45a0': 'F',
        'd4c2': 'O',
        'fb1c': 'H',
        '529f': 'I',
        '8f61': 'G',
        '0fa7': 'T',
        '09be': 'N',
        '05d9': 'S',
        '644f': 'C',
        'bdb9': 'L',
        '237f': 'D',
        '7ba1': 'U',
        '591e': 'I', # Assuming this is a valid substitution despite 'I' already being mapped.
        'afc5': 'A',
        '3a8d': 'C',
    }

    # The encrypted words
    words_hex = [
        "45a0afc53a8dafc50fa7529f",
        "7ba1d4c2644ffb1c05d9",
        "bdb9d4c2591e05d9529f05d9",
        "d4c245a0",
        "644fd4c2bdb9237f",
        "fb1c529f8f61fb1c",
        "fb1c0fa7529f8f6109be05d9"
    ]

    deciphered_phrase = []
    for hex_word in words_hex:
        deciphered_word = ""
        # Process each word in 4-symbol (8-char) chunks
        for i in range(0, len(hex_word), 8):
            chunk = hex_word[i:i+8]
            # Use a slice of the hex code since the original problem statement used 4 symbols, not 8 hex chars
            # But since all are unique, we can use the 8-char hex string as the key.
            # Let's adjust keys to be 4-char for clarity with the problem.
            # No, the hex codes are the 'numbers'. Let's stick to the full code.
            key = hex_word[i:i+4]
            # It seems the original problem had a typo and used variable length hex codes.
            # Let's re-read "several symbols mean one number". Hex values are `45a0`, `afc5`, etc.
            # My assumption of 4-symbol chunks was based on 8 hex chars / 2 letters. This is 4 hex chars.
            # Let's adjust the code to use 4-char chunks.
            chunk_key = hex_word[i:i+4]
            if chunk_key in mapping:
                 deciphered_word += mapping[chunk_key]
            else:
                 deciphered_word += "?" # Placeholder for any unmapped code
        deciphered_phrase.append(deciphered_word)
        
    # Let's rebuild the mapping and phrase based on the correct 4-char chunking
    # This was a mistake in my python implementation plan, the logic remains sound.
    
    final_mapping = {
        '45a0': 'F', 'afc5': 'A', '3a8d': 'C', '0fa7': 'T', '529f': 'I',
        '7ba1': 'U', 'd4c2': 'O', '644f': 'C', 'fb1c': 'H', '05d9': 'S',
        'bdb9': 'L', '591e': 'U',
        '237f': 'D',
        '8f61': 'G',
        '09be': 'N'
    }

    final_phrase_list = []
    all_words_hex = [
        "45a0afc53a8dafc50fa7529f", "7ba1d4c2644ffb1c05d9", "bdb9d4c2591e05d9529f05d9",
        "d4c245a0", "644fd4c2bdb9237f", "fb1c529f8f61fb1c", "fb1c0fa7529f8f6109be05d9"
    ]

    for word in all_words_hex:
        decoded_word = ""
        for i in range(0, len(word), 4):
            key = word[i:i+4]
            decoded_word += final_mapping.get(key, "?")
        final_phrase_list.append(decoded_word)

    final_phrase = " ".join(final_phrase_list)

    print("Deciphered Phrase:")
    print(final_phrase)
    print("\nSubstitution Map (Equation):")
    # Print each substitution on a new line
    for code, letter in sorted(final_mapping.items()):
        print(f"{code} = {letter}")

solve_crypto()