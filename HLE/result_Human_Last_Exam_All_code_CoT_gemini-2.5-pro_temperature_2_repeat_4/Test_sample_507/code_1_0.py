def solve_linguistic_puzzle():
    """
    Identifies and prints the two Asian languages with similar words for "mom", "dad", and "broom".
    """
    lang1_name = "Korean"
    lang1_culture = "Korean (East Asia)"
    lang2_name = "Tamil"
    lang2_culture = "Tamil (South Asia, Dravidian)"

    word1 = "mom"
    lang1_word1 = "eomma (엄마)"
    lang2_word1 = "amma (அம்மா)"

    word2 = "dad"
    lang1_word2 = "appa (아빠)"
    lang2_word2 = "appa (அப்பா)"

    word3 = "broom"
    lang1_word3 = "bitjaru (빗자루)"
    lang2_word3 = "vilakkumaaru (விளக்குமாறு)"

    print("This puzzle requires finding two languages from different Asian cultures with no direct contact that share similar words.")
    print("The chosen languages are Korean and Tamil.\n")
    print(f"Language 1: {lang1_name} (Culture: {lang1_culture})")
    print(f"Language 2: {lang2_name} (Culture: {lang2_culture})")
    print("-" * 50)
    print("Here is the comparison of the words:\n")

    # --- MOM ---
    print(f"WORD: '{word1.upper()}'")
    print(f"{lang1_name}: {lang1_word1}")
    print(f"{lang2_name}:   {lang2_word1}")
    print("Similarity Analysis: The words are nearly identical, both using the simple labial 'm' sound typical of nursery words for 'mother'.\n")

    # --- DAD ---
    print(f"WORD: '{word2.upper()}'")
    print(f"{lang1_name}: {lang1_word2}")
    print(f"{lang2_name}:   {lang2_word2}")
    print("Similarity Analysis: The words are identical. This is a well-known linguistic coincidence.\n")

    # --- BROOM ---
    print(f"WORD: '{word3.upper()}'")
    print(f"{lang1_name}: {lang1_word3}")
    print(f"{lang2_name}:   {lang2_word3}")
    print("Similarity Analysis: While phonetically different overall, a notable similarity exists.")
    print("Both languages form a compound word for 'broom', and the second parts, '-jaru' (meaning handle/tool) and '-maaru' (a sweeping tool), are phonetically similar, especially their 'aru' ending.")

if __name__ == '__main__':
    solve_linguistic_puzzle()
<<<Korean and Tamil>>>