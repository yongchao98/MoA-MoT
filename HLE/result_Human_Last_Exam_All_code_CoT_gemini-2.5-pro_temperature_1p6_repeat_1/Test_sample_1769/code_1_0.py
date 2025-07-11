def find_native_word_for_sister():
    """
    This script calculates what the word for "sister" would be in Modern English
    if it had evolved from its Old English root without influence from Old Norse.
    """

    # Step 1: State the current word and its origin.
    current_word = "sister"
    origin_language = "Old Norse"
    origin_word = "systir"
    print(f"The current Modern English word '{current_word}' is not native to English.")
    print(f"It was adopted from the {origin_language} word '{origin_word}' during the Viking invasions (c. 8th-11th centuries).\n")

    # Step 2: Identify the original Old English word.
    old_english_word = "sweostor"
    print(f"The original Anglo-Saxon (Old English) word was '{old_english_word}'.\n")

    # Step 3: Explain the hypothetical evolution.
    print("Without the Norse influence, 'sweostor' would have undergone natural sound changes:")
    
    middle_english_form = "suster"
    print(f" - First, 'sweostor' simplified in Middle English to a form like '{middle_english_form}'.")
    print("   (This is similar to how German has 'Schwester' and Dutch has 'zuster').")

    final_word = "suster"
    print(f" - This Middle English form '{middle_english_form}' would then evolve into the Modern English word '{final_word}'.")
    print("   (The 'u' vowel sound would be pronounced like the 'u' in 'cup' or 'sun').\n")

    # Step 4: Final Conclusion
    print("Therefore, if the Norse had never invaded, the Modern English word for a female sibling would be:")
    print(final_word)

if __name__ == '__main__':
    find_native_word_for_sister()
    print("\n<<<suster>>>")