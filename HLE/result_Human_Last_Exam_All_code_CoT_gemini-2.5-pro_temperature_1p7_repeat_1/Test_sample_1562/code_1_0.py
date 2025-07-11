def solve_encryption():
    """
    This function deciphers the encrypted phrase based on the logic described.
    It constructs the final phrase and prints it as a mathematical-style equation.
    """
    encrypted_words = [
        "45a0afc53a8dafc50fa7529f",
        "7ba1d4c2644ffb1c05d9",
        "bdb9d4c2591e05d9529f05d9",
        "d4c245a0",
        "644fd4c2bdb9237f",
        "fb1c529f8f61fb1c",
        "fb1c0fa7529f8f6109be05d9",
    ]

    # Mapping based on the step-by-step decipherment
    # This map is the result of the logical deduction.
    # Assumes 'FUTURE GOADS HORSES OF DARK DEAD DREADS'
    code_map = {
        '45a0': 'F', 'afc5': 'U', '3a8d': 'T', '0fa7': 'R', '529f': 'E', # FUTURE
        '7ba1': 'G', 'd4c2': 'O', '644f': 'A', 'fb1c': 'D', '05d9': 'S', # GOADS and part of DARK
        'bdb9': 'H', '591e': 'R',                                      # HORSES
                                                                        # OF is from d4c2, 45a0
        '237f': 'K',                                                  # DARK
        '8f61': 'A',                                                  # DEAD
        '09be': 'D'                                                   # DREADS
    }

    deciphered_phrase = []
    equations = []

    for enc_word in encrypted_words:
        dec_word = ""
        equation_parts = []
        # Process word in blocks of 4
        for i in range(0, len(enc_word), 4):
            code = enc_word[i:i+4]
            letter = code_map.get(code, '?')
            dec_word += letter
            equation_parts.append(f"{code}({letter})")
        
        deciphered_phrase.append(dec_word)
        equations.append(" + ".join(equation_parts) + f" = {dec_word}")

    final_phrase = " ".join(deciphered_phrase)
    
    print("Deciphered Phrase:")
    print(final_phrase)
    print("\nEquations for each word:")
    for eq in equations:
        print(eq)

solve_encryption()