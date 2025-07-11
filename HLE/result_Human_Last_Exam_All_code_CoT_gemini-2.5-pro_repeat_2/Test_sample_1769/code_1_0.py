def find_hypothetical_word():
    """
    This function explains and determines what the Modern English word for 'sister'
    would be if the Norse had not influenced the English language.
    """

    # The current word and its origin
    current_word = "sister"
    origin = "Old Norse 'systir'"

    # The original native word
    old_english_word = "sweoster"

    # Cognates in related languages
    dutch_cognate = "zuster"
    german_cognate = "Schwester"

    # The likely evolution path
    # A common sound change in English was 'swe-' becoming 'su-'.
    # For example, Old English 'swilc' became Modern English 'such'.
    # Applying this pattern to 'sweoster' gives 'suster'.
    hypothetical_word = "suster"

    print("--- The Etymological Path of 'Sister' ---")
    print(f"1. The modern word '{current_word}' was borrowed into English from {origin}.")
    print(f"2. Before the Norse influence, the native Old English (Anglo-Saxon) word was '{old_english_word}'.")
    print(f"3. Without the borrowing, '{old_english_word}' would have evolved naturally.")
    print(f"4. A common linguistic development in English was the change from 'swe-' to 'su-', as seen in Old English 'swilc' becoming 'such'.")
    print(f"5. Applying this same change to '{old_english_word}' would have resulted in '{hypothetical_word}'.")
    print(f"6. This is supported by similar words in other West Germanic languages, like Dutch '{dutch_cognate}'.")
    print("\nTherefore, the likely Modern English word for 'sister' would be:")
    print(f"'{hypothetical_word}'")

find_hypothetical_word()