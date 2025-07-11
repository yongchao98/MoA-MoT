def find_similar_words():
    """
    This function presents two languages from distinct Asian cultures
    that have remarkably similar words for "mom," "dad," and "broom."
    """
    language1 = "Korean"
    culture1 = "East Asia (Koreanic language family)"

    language2 = "Tamil"
    culture2 = "South Asia (Dravidian language family)"

    print("The two languages are Korean and Tamil.")
    print("These languages are from two distinct Asian cultures with no direct genealogical relationship or significant historical contact.")
    print("-" * 30)

    # Word for Mom
    mom_ko = "엄마 (eomma)"
    mom_ta = "அம்மா (ammā)"
    print(f"Word for 'Mom':")
    print(f"  In {language1}: {mom_ko}")
    print(f"  In {language2}: {mom_ta}")
    print("  Similarity: Nearly identical. Both are based on a simple vowel + 'mm' sound.\n")

    # Word for Dad
    dad_ko = "아빠 (appa)"
    dad_ta = "அப்பா (appā)"
    print(f"Word for 'Dad':")
    print(f"  In {language1}: {dad_ko}")
    print(f"  In {language2}: {dad_ta}")
    print("  Similarity: Nearly identical. Both are based on a simple vowel + 'pp' sound.\n")

    # Word for Broom
    broom_ko_full = "빗자루 (bitjaru)"
    broom_ko_root = "비 (bi)"
    broom_ta = "விளக்குமாறு (viḷakkumāṟu)"
    print(f"Word for 'Broom':")
    print(f"  In {language1}: {broom_ko_full}, where the root for 'broom' is '{broom_ko_root}'.")
    print(f"  In {language2}: {broom_ta}")
    print("  Similarity: This is more subtle. The Korean root 'bi' and the initial sound of the Tamil word 'vi-' are phonetically very close (a voiced bilabial plosive vs. a voiced labio-dental approximant), a similarity found in a surprising number of coincidences between unrelated languages.")

find_similar_words()