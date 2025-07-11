def solve_riddle():
    """
    This function solves the riddle by connecting Vienna's boulevards
    to an English poet's name through a Russian translation.
    """

    # Step 1: Identify the key location in Vienna.
    vienna_boulevard_german = "Ringstrasse"
    key_german_word = "Ring"
    
    print(f"The wide boulevards of Vienna are known as the '{vienna_boulevard_german}'.")
    print(f"The key German word in the name is '{key_german_word}'.")

    # Step 2: Find the Russian translation.
    russian_word_for_ring = "кольцо"
    pronunciation = "kol'tso"
    
    print(f"The Russian word for 'ring' is '{russian_word_for_ring}', pronounced roughly as '{pronunciation}'.")

    # Step 3: Find the English poet with a similar-sounding surname.
    english_poet_surname = "Coleridge"
    russian_transliteration = "Кольридж"
    
    print(f"The English poet's surname '{english_poet_surname}' is written as '{russian_transliteration}' in Russian.")
    print(f"This creates a phonetic pun with the Russian word for ring '{russian_word_for_ring}'.")
    
    # Final Answer
    print("\nTherefore, the surname of the English poet is:")
    print(english_poet_surname)

solve_riddle()