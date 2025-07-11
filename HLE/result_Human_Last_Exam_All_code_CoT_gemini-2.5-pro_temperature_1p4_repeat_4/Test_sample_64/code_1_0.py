def find_pitch_accent():
    """
    This function looks up and explains the pitch accent for the word 「弟」.
    """
    word = "弟"
    reading = "おとうと"
    morae = ["お", "と", "う", "と"]

    # In standard Japanese pitch accent dictionaries, 「弟」 is listed with
    # the pattern [0].
    pitch_accent_number = 0
    
    # The [0] pattern is called Heiban (平板).
    pattern_name = "Heiban"
    
    print(f"Analyzing the pitch accent for: 「{word}」 ({reading})")
    print("-" * 30)
    print(f"The word is read as '{reading}', which consists of {len(morae)} morae: {', '.join(morae)}.")
    
    # Heiban [0] pattern explained:
    # The pitch starts low on the first mora, rises on the second, and
    # remains high for the rest of the word and any following particles.
    print(f"\nThe dictionary pitch accent is [{pitch_accent_number}], which is known as the '{pattern_name}' pattern.")
    print("\nThe pitch contour is as follows:")
    print("Mora 1 (お): Low")
    print("Mora 2 (と): High")
    print("Mora 3 (う): High")
    print("Mora 4 (と): High")
    print("\nThis pattern corresponds to choice A.")

if __name__ == "__main__":
    find_pitch_accent()