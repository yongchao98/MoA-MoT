def solve_riddle():
    """
    This function explains the steps to solve the riddle about the
    English poet mentioned in a Russian text about Vienna.
    """
    # Step 1: Identify the famous boulevard in Vienna.
    vienna_boulevard = "Ringstrasse"
    print(f"Step 1: The famous wide boulevard in Vienna is the '{vienna_boulevard}'.")

    # Step 2: Translate the key part of the name into Russian.
    german_word = "Ring"
    russian_equivalent_adjective = "Кольцевой"
    russian_pronunciation = "Koltsevoy"
    print(f"Step 2: The German word '{german_word}' in '{vienna_boulevard}' corresponds to the Russian adjective '{russian_equivalent_adjective}' (pronounced roughly as '{russian_pronunciation}').")

    # Step 3: Find the English poet whose name creates a pun in Russian.
    poet_surname_english = "Coleridge"
    poet_surname_russian = "Кольридж"
    poet_pronunciation = "Kol'ridzh"
    print(f"Step 3: The surname of the English poet '{poet_surname_english}' is transliterated into Russian as '{poet_surname_russian}' (pronounced roughly as '{poet_pronunciation}').")

    # Step 4: Explain the pun and provide the final answer.
    print(f"Step 4: The pun is based on the phonetic similarity between the word for the 'Ring' boulevard ('{russian_equivalent_adjective}') and the poet's surname ('{poet_surname_russian}').")
    
    print("\nTherefore, the surname is:")
    print(poet_surname_english)

solve_riddle()