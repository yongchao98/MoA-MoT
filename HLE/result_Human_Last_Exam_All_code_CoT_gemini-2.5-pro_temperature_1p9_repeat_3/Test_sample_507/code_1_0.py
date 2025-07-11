def find_similar_words():
    """
    Identifies and prints a comparison of two Asian languages from different
    cultures with no direct contact, where words for "mom", "dad", and "broom"
    show surprising similarities.
    """

    lang1 = "Korean"
    culture1 = "East Asia"

    lang2 = "Telugu"
    culture2 = "South-Central India (Dravidian)"

    print(f"Finding similar words in two languages from distant Asian cultures:\n")
    print(f"Language 1: {lang1} ({culture1})")
    print(f"Language 2: {lang2} ({culture2})")
    print("-" * 50)
    print("These two languages have no direct historical contact and belong to different language families (Koreanic and Dravidian).\n")

    print("WORD COMPARISON:\n")

    # --- Mom ---
    print("1. Word: MOM")
    print(f"   In {lang1}: 엄마 (eomma)")
    print(f"   In {lang2}: అమ్మ (amma)")
    print("   Analysis: The words are nearly identical.\n")

    # --- Dad ---
    print("2. Word: DAD")
    print(f"   In {lang1}: 아빠 (appa)")
    print(f"   In {lang2}: అప్ప (appa - a word for father)")
    print("   Analysis: The words are identical.\n")

    # --- Broom ---
    print("3. Word: BROOM")
    print(f"   In {lang1}: 빗자루 (bitjaru)")
    print(f"   In {lang2}: చీపురు (cheepuru)")
    print("   Analysis: While not identical, these words share a similar structure and phonetic elements. Both are trisyllabic and share the same final '-ru' sound, with the initial consonant sounds /b/ and /ch/ also being phonetically related (plosive vs. affricate). This represents a notable, if less exact, similarity.")

if __name__ == '__main__':
    find_similar_words()