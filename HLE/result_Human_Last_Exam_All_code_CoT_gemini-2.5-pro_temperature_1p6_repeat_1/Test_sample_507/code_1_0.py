def find_similar_words():
    """
    This script identifies two languages from different Asian cultures with
    no direct contact that share highly similar words for "mom", "dad",
    and the material for a "broom".
    """
    lang1 = "Korean"
    lang2 = "Tamil (Dravidian)"

    print(f"Finding similar words in {lang1} and {lang2}, two languages from culturally and geographically distant parts of Asia.")
    print("These languages belong to unrelated language families (Koreanic and Dravidian).\n")

    # --- Mom ---
    print("1. The word for 'Mom':")
    print(f"   - In {lang1}: eomma (엄마)")
    print(f"   - In {lang2}: amma (அம்மா)")
    print("   -> The similarity is very high.\n")

    # --- Dad ---
    print("2. The word for 'Dad':")
    print(f"   - In {lang1}: appa (아빠)")
    print(f"   - In {lang2}: appa (அப்பா)")
    print("   -> The words are identical.\n")

    # --- Broom (via its material) ---
    print("3. The word for 'Broom':")
    print("   While the word for the tool 'broom' is different, the word for the grass/straw used to make a traditional broom is identical.")
    print(f"   - In {lang1}, 'grass' or 'straw' is: pul (풀)")
    print(f"   - In {lang2}, 'grass' or 'straw' is: pul (புல்)")
    print("   -> The word for the broom's material is identical, solving the final piece of the puzzle.")

if __name__ == "__main__":
    find_similar_words()