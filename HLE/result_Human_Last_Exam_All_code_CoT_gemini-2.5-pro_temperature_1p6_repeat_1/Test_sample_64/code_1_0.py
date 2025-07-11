def find_pitch_accent():
    """
    This function demonstrates the pitch accent for the Japanese word 「弟」.
    It uses a dictionary as a simple knowledge base.
    """
    # A dictionary mapping words to their pitch accent information.
    # The format is: Kanji: (Reading, Pitch Accent Number, Pattern Name)
    # Pitch Accent Numbers:
    # 0: Heiban (flat)
    # 1: Atamadaka (high on first mora, then drops)
    # >1: Nakadaka (high on a middle mora) or Odaka (high until the last mora, then drops)
    pitch_accent_db = {
        "弟": ("おとうと", 0, "Heiban"),
        "頭": ("あたま", 2, "Nakadaka"),
        "犬": ("いぬ", 2, "Odaka"), # In this case it's Odaka because it's the last mora
        "箸": ("はし", 1, "Atamadaka"),
        "橋": ("はし", 2, "Odaka")
    }

    word_to_find = "弟"
    
    if word_to_find in pitch_accent_db:
        reading, pattern_number, pattern_name = pitch_accent_db[word_to_find]
        
        print(f"The Japanese word is 「{word_to_find}」, which is read as 「{reading}」.")
        print("The standard pitch accent pattern can be represented by a number.")
        print(f"The equation for this word's pitch accent is:")
        print(f"「{word_to_find}」 = [{pattern_number}]")
        print(f"\nA pitch accent of [{pattern_number}] corresponds to the '{pattern_name}' pattern.")
        print("For a Heiban word like 「弟」 (おとうと), the pitch is low on the first mora (お) and then rises and stays high for the rest of the word (とうと).")

    else:
        print(f"Pitch accent information for '{word_to_find}' not found.")

find_pitch_accent()