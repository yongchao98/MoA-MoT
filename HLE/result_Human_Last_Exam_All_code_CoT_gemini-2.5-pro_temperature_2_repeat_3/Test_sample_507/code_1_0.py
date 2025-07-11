def find_similar_words():
    """
    Prints a comparison of words between Korean and Tamil, two unrelated
    Asian languages with surprising lexical similarities.
    """
    language_comparison = {
        "Word": ("Korean", "Tamil"),
        "Mom": ("엄마 (eomma)", "அம்மா (ammā)"),
        "Dad": ("아빠 (appa)", "அப்பா (appā)"),
        "Broom": ("빗자루 (bitjaru)", "துடைப்பம் (thudaippam)")
    }

    print("Comparing words in Korean and Tamil:")
    print("-" * 50)
    
    # Header
    word, (lang1, lang2) = language_comparison.pop("Word")
    print(f"{'Word':<10} | {lang1:<20} | {lang2:<20}")
    print("-" * 50)

    # Print each word pair
    for word_type, (word1, word2) in language_comparison.items():
        print(f"{word_type:<10} | {word1:<20} | {word2:<20}")
    
    print("-" * 50)
    print("\nObservation:")
    print("- 'Mom' and 'Dad' are nearly identical in both languages.")
    print("- 'Broom' does not show a direct similarity.")
    print("- This illustrates a remarkable linguistic coincidence between two distant, unrelated cultures.")

# Execute the function to display the results
find_similar_words()