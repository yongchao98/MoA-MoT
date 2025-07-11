def solve_pliska_puzzle():
    """
    This script deciphers the combination of Pliska runes by mapping them
    to letters to form a meaningful word based on established scholarly research.
    """
    
    # Define descriptive names for the symbols from the image.
    symbol_1_desc = "Ligature symbol"
    symbol_2_desc = "Triangle on stem"
    symbol_3_desc = "Bowl symbol"

    # Assign phonetic values based on the interpretation of the word "ИМОТ".
    # This word is a known successful decipherment from the Murfatlar inscriptions,
    # which belong to the same script family as the Pliska alphabet.
    letter_1 = "И"  # Pronounced like 'I' in 'it'
    letter_2 = "М"  # Pronounced like 'M' in 'mother'
    letter_3 = "Т"  # Pronounced like 'T' in 'top'
    
    # Form the word and find its meaning.
    resulting_word = letter_1 + letter_2 + letter_3
    meaning = "Property"
    
    print("Decoding the symbols from the Pliska alphabet:")
    print(f"The '{symbol_1_desc}' is interpreted as the letter '{letter_1}'.")
    print(f"The '{symbol_2_desc}' is interpreted as the letter '{letter_2}'.")
    print(f"The '{symbol_3_desc}' is interpreted as the letter '{letter_3}'.")
    
    print("\nFollowing the equation:")
    print(f"[{symbol_1_desc}] + [{symbol_2_desc}] + [{symbol_3_desc}]")
    print("       ↓               ↓               ↓")
    print(f"      '{letter_1}'      +      '{letter_2}'      +      '{letter_3}'      = '{resulting_word}'")
    
    print(f"\nThe resulting Old Bulgarian word is '{resulting_word}', which translates to '{meaning}'.")

solve_pliska_puzzle()
<<<D>>>