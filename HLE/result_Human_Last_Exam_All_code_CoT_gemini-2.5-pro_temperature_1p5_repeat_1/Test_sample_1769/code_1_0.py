def find_native_word_for_sister():
    """
    This function explains and prints the likely Modern English word for 'sister'
    if the Norse language had not influenced English.
    """
    
    # 1. Identify the original Old English and competing Old Norse words.
    old_english_word = "sweostor"
    old_norse_word = "systir"

    print("--- The Linguistic Problem ---")
    print(f"Our current word 'sister' comes from the Old Norse word '{old_norse_word}'.")
    print(f"To find the word without Norse influence, we must trace the evolution of the native Old English word: '{old_english_word}'.")
    print("-" * 30)

    # 2. Explain the evolution from Old English to Middle English.
    middle_english_form = "suster"
    
    print("\n--- The Likely Evolution ---")
    print(f"The Old English '{old_english_word}' evolved naturally into the Middle English form '{middle_english_form}'.")
    print("This happened through common sound changes:")
    print("1. Simplification of the initial 'sw-' sound to 's-'.")
    print("2. A shift in the vowel sound 'eo' to 'u'.")
    print("-" * 30)

    # 3. Present the conclusion.
    final_word = middle_english_form
    
    print("\n--- The Conclusion ---")
    print(f"The word '{final_word}' was widely used in parts of England for centuries.")
    print("Without the competition from the Norse '{old_norse_word}', it is the most likely candidate for the modern word.")
    print("\nPredicted Modern English word for 'sister':")
    print(final_word)

find_native_word_for_sister()