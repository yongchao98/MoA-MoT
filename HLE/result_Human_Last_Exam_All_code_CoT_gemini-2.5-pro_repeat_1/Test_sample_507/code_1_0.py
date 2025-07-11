def find_similar_words():
    """
    This function presents a solution to the linguistic puzzle by identifying
    two languages from different Asian cultures with similar words for "mom",
    "dad", and "broom".
    """
    lang1 = "Korean"
    lang2 = "Persian (Farsi)"
    
    print(f"The two languages are {lang1} and {lang2}.")
    print("These languages are from entirely different language families and have no history of direct contact, making the similarities a remarkable coincidence.\n")

    # --- Mom ---
    mom_lang1 = "eomma (엄마)"
    mom_lang2 = "mâmân (مامان)"
    print("1. Word for 'Mom':")
    print(f"   - In {lang1}: {mom_lang1}")
    print(f"   - In {lang2}: {mom_lang2}")
    print("   - Similarity: Both are two-syllable words starting with the 'm' sound, a common feature for 'mom' words globally (nursery words).\n")

    # --- Dad ---
    dad_lang1 = "appa (아빠)"
    dad_lang2 = "bâbâ (بابا)"
    print("2. Word for 'Dad':")
    print(f"   - In {lang1}: {dad_lang1}")
    print(f"   - In {lang2}: {dad_lang2}")
    print("   - Similarity: Both words are formed by reduplicating a simple syllable with a labial consonant ('p' or 'b') and the 'a' vowel, another classic nursery word pattern.\n")

    # --- Broom ---
    broom_lang1 = "bitjaru (빗자루)"
    broom_lang2 = "jâru (جارو)"
    print("3. Word for 'Broom':")
    print(f"   - In {lang1}: {broom_lang1}")
    print(f"   - In {lang2}: {broom_lang2}")
    print("   - Similarity: This is the most striking and unlikely coincidence. The Korean word is a compound, but its second part, 'jaru' (자루), is phonetically almost identical to the complete Persian word 'jâru'.")

if __name__ == "__main__":
    find_similar_words()
