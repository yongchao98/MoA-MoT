def find_similar_language_words():
    """
    Identifies and prints information about two Asian languages with no direct
    contact that share similar words for "mom", "dad", and "broom".
    """
    language1 = "Tagalog"
    culture1 = "Philippines (Austronesian Language Family)"
    language2 = "Tamil"
    culture2 = "Southern India / Sri Lanka (Dravidian Language Family)"

    print(f"The two languages are {language1} and {Tamil}.")
    print("They come from completely different language families and geographically distant cultures with no modern direct contact.\n")

    print("--- Word Comparison ---\n")

    # --- Mom ---
    mom1_word = "Nanay"
    mom2_word = "Annai (அன்னை)"
    print(f"1. The word for 'Mom':")
    print(f"   - In {language1}: {mom1_word}")
    print(f"   - In {language2}: {mom2_word}")
    print("   Note: 'Annai' is a formal/poetic term for mother in Tamil, making the similarity to 'Nanay' very striking.\n")

    # --- Dad ---
    dad1_word = "Tatay"
    dad2_word = "Thathai (தாதை)"
    print(f"2. The word for 'Dad':")
    print(f"   - In {language1}: {dad1_word}")
    print(f"   - In {language2}: {dad2_word}")
    print("   Note: 'Thathai' is a formal/poetic term for father in Tamil, which is phonetically very similar to 'Tatay'.\n")

    # --- Broom ---
    broom1_word = "Walis"
    broom2_word = "Vāri (வாரி)"
    print(f"3. The word for 'Broom':")
    print(f"   - In {language1}: {broom1_word}")
    print(f"   - In {language2}, the related verb 'to sweep/gather up' is: {broom2_word}")
    print("   Note: The Tagalog 'Walis' refers to both the broom and the act of sweeping. It is phonetically similar to the Tamil verb 'Vāri', which means to sweep or scrape together.\n")

if __name__ == '__main__':
    find_similar_language_words()