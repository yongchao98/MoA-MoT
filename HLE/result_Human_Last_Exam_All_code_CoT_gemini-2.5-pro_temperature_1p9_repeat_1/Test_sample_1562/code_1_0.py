def solve_cipher():
    """
    This function deciphers the encrypted phrase based on a derived mapping
    and prints the final phrase with the corresponding numbers for each letter.
    """
    
    # The derived mapping from symbol chunks to letters
    mapping = {
        '45a0': 'F', 'afc5': 'A', '3a8d': 'I', '0fa7': 'L', '529f': 'E',
        '7ba1': 'T', 'd4c2': 'O', '644f': 'U', 'fb1c': 'C', '05d9': 'S',
        'bdb9': 'H', '591e': 'R', '237f': 'Y', '8f61': 'K', '09be': 'G'
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
    
    final_phrase_parts = []
    
    for word_hex in words_hex:
        deciphered_word = ""
        equation_parts = []
        # Process each word in 4-character chunks
        for i in range(0, len(word_hex), 4):
            chunk = word_hex[i:i+4]
            letter = mapping.get(chunk, '?')
            deciphered_word += letter
            equation_parts.append(f"{chunk}({letter})")
        final_phrase_parts.append(f"{deciphered_word} [{ ' + '.join(equation_parts) }]")

    print("\n".join(final_phrase_parts))

solve_cipher()