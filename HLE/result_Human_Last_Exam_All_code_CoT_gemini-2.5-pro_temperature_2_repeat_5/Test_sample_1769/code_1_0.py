def find_hypothetical_word():
    """
    Deduces the modern English word for 'sister' without Norse influence
    by tracing the etymology of the original Old English word.
    """
    old_english_word = "sweostor"
    middle_english_form = "suster"
    
    print("This script determines the Modern English word for 'sister' had the Norse never introduced their term.")
    print("-" * 20)
    
    # Step 1: Explain the origin of the current word
    print("1. Our current word, 'sister', comes from the Old Norse word 'systir'.")
    print("-" * 20)

    # Step 2: Identify the native word that would have been used instead
    print(f"2. Without that influence, the word would have evolved from the Old English (Anglo-Saxon) word: '{old_english_word}'.")
    print("-" * 20)

    # Step 3: Show the evolution path
    print(f"3. Let's trace the hypothetical evolution of '{old_english_word}':")
    print(f"   - From Old to Middle English, '{old_english_word}' commonly became '{middle_english_form}'.")
    print("     This is due to a known sound change where 'sw-' at the start of a word simplified to 's-' (e.g., Old English 'swa' became 'so').")
    print("-" * 20)
    
    # Step 4: Project the Modern English result
    print(f"4. The Middle English form '{middle_english_form}' would then have evolved into a Modern English word.")
    print("     - The pronunciation would have changed during the Great Vowel Shift, similar to 'butter' or 'utter'.")
    print("     - The spelling would have likely been standardized to 'suster'.")
    print("-" * 20)

    # Step 5: State the final conclusion
    final_word = middle_english_form
    print("The final equation is essentially: ")
    print(f"{old_english_word} (Old English) -> {middle_english_form} (Middle English) -> {final_word} (Modern English)")
    print("\nTherefore, the likely Modern English word for 'sister' would be:")
    print(final_word)

find_hypothetical_word()