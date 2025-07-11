def find_hypothetical_word():
    """
    This script reconstructs the Modern English word for "sister"
    by tracing its Old English root through historical sound changes,
    ignoring the actual borrowing from Old Norse.
    """
    print("--- Linguistic Reconstruction: The Word for 'Sister' ---")

    # Step 1: State the origin of the current word.
    current_word = "sister"
    origin = "Old Norse 'systir'"
    print(f"\n1. The current word '{current_word}' is not native to English.")
    print(f"   It was borrowed from {origin}, replacing the original English word.")

    # Step 2: Identify the original Old English word.
    old_english_word = "sweostor"
    print(f"\n2. The native Old English word was '{old_english_word}'.")

    # Step 3: Trace the evolution from Old to Middle English.
    # The 'eo' diphthong in 'sweostor' commonly evolved into 'u'.
    # The unstressed final syllable '-or' weakened to '-er'.
    # This process gives us the intermediate Middle English form.
    middle_english_word = "suster"
    print("\n3. Evolution from Old English to Middle English:")
    print(f"   Old English '{old_english_word}' -> Middle English '{middle_english_word}'")
    print("   - The 'eo' sound shifted to 'u'.")
    print("   - The unstressed ending '-or' became '-er'.")
    print("   (Note: The form 'suster' was widely used in Middle English).")

    # Step 4: Trace the evolution from Middle to Modern English.
    # The primary change here is the "short u" vowel sound /ʊ/ shifting to /ʌ/,
    # as seen in the evolution of 'sunne' to 'sun' or 'cutte' to 'cut'.
    modern_hypothetical_word = "suster"
    print("\n4. Evolution from Middle English to Modern English:")
    print(f"   Middle English '{middle_english_word}' -> Modern English '{modern_hypothetical_word}'")
    print("   - The vowel 'u' (pronounced like in 'put') shifted to the sound in 'cup'.")
    print("   - The spelling would likely have remained the same.")

    # Step 5: Final Conclusion
    print("\n--- Conclusion ---")
    print("Without the Norse influence, the Old English 'sweostor' would have naturally evolved through the form 'suster'.")
    print(f"Therefore, the Modern English word for 'sister' would most likely be:")
    print(f"\n'{modern_hypothetical_word}'")

find_hypothetical_word()