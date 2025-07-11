def decode_pliska_symbols():
    """
    Decodes a combination of Pliska alphabet symbols based on Vasil Ĭonchev's research.
    """
    # According to Vasil Ĭonchev's interpretation of Proto-Bulgarian runes:
    # The first symbol corresponds to the Cyrillic letter 'Д' (D).
    # The second symbol corresponds to the Cyrillic letter 'А' (A).
    # The third symbol corresponds to the Cyrillic letter 'Р' (R).
    
    symbol_letters = {
        "first_symbol": "Д",
        "second_symbol": "А",
        "third_symbol": "Р"
    }

    letter1 = symbol_letters["first_symbol"]
    letter2 = symbol_letters["second_symbol"]
    letter3 = symbol_letters["third_symbol"]

    # Combine the letters to form the word
    word_in_cyrillic = letter1 + letter2 + letter3
    
    # Translate the word to English
    english_translation = "Gift"

    print("Decoding the symbols from the Pliska alphabet according to Vasil Ĭonchev's research:")
    print(f"The first symbol represents the letter: {letter1}")
    print(f"The second symbol represents the letter: {letter2}")
    print(f"The third symbol represents the letter: {letter3}")
    print("\nForming the word from the letters in the equation:")
    print(f"'{letter1}' + '{letter2}' + '{letter3}' => '{word_in_cyrillic}'")
    print(f"\nThe Bulgarian word '{word_in_cyrillic}' translates to '{english_translation}' in English.")

decode_pliska_symbols()