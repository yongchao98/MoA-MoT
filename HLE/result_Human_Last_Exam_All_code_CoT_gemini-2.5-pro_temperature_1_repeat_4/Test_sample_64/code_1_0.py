def get_pitch_accent_info(word, hiragana, morae, accent_pattern_numerical):
    """
    Determines and prints the pitch accent pattern and its classification.
    accent_pattern_numerical: 0 for Heiban, 1 for Atamadaka, >1 for Nakadaka/Odaka
                              based on where the drop occurs. For おとうと, it's 3.
    """
    num_morae = len(morae)
    pitch_pattern = []

    # Generate L/H pattern from numerical representation
    if accent_pattern_numerical == 0:  # Heiban
        pattern_name = "Heiban"
        pitch_pattern.append("L")
        for _ in range(num_morae - 1):
            pitch_pattern.append("H")
    elif accent_pattern_numerical == 1:  # Atamadaka
        pattern_name = "Atamadaka"
        pitch_pattern.append("H")
        for _ in range(num_morae - 1):
            pitch_pattern.append("L")
    else:  # Nakadaka or Odaka
        if accent_pattern_numerical < num_morae + 1: # Nakadaka
             pattern_name = "Nakadaka"
        else: # This case is not for this word.
             pattern_name = "Odaka"

        pitch_pattern.append("L") # First mora is always low (unless Atamadaka)
        for i in range(1, num_morae):
            if i < accent_pattern_numerical:
                pitch_pattern.append("H")
            else:
                pitch_pattern.append("L")

    # --- Output ---
    print(f"Word: 「{word}」 ({hiragana})")
    print("Pitch Accent Analysis:")
    for i in range(num_morae):
        print(f"  Mora {i+1} 「{morae[i]}」: {pitch_pattern[i]}")

    print(f"\nThis pattern, where the pitch drops within the word, is classified as: {pattern_name}")
    print("\nAnswer Choice C corresponds to Nakadaka.")


# Standard data for the word 「弟」
word_kanji = "弟"
word_hiragana = "おとうと"
word_morae = ["お", "と", "う", "と"]
# The accent is [3], meaning the pitch drops after the 3rd mora.
accent_core_position = 3

get_pitch_accent_info(word_kanji, word_hiragana, word_morae, accent_core_position)