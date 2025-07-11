def print_pitch_accent_info(kanji, kana, pattern_name, pattern_code):
    """
    Prints information about a Japanese word's pitch accent.

    Args:
        kanji (str): The word in Kanji.
        kana (str): The word in Hiragana.
        pattern_name (str): The name of the pitch accent pattern (e.g., "Heiban").
        pattern_code (str): A simple representation of the pitch (e.g., "L-H-H-H").
    """
    print(f"Word: 「{kanji}」 ({kana})")
    print(f"Standard Pitch Accent Pattern: {pattern_name}")
    print("\nThis means the pitch contour is as follows:")

    moras = ["お", "と", "う", "と"]
    pitches = pattern_code.split('-')

    for i in range(len(moras)):
        print(f"Mora {i+1}「{moras[i]}」: {pitches[i]}")

    print("\nWhen followed by a particle like が, the particle also stays high:")
    print("お(L)とう(H)と(H)が(H)")


# Information for the word 「弟」
word_kanji = "弟"
word_kana = "おとうと"
accent_pattern = "Heiban (平板)"
pitch_representation = "L-H-H-H" # L=Low, H=High

# Execute the function to display the answer
print_pitch_accent_info(word_kanji, word_kana, accent_pattern, pitch_representation)