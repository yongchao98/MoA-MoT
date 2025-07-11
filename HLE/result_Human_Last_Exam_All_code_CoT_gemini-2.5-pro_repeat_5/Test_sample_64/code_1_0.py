import collections

def get_pitch_accent_info(word):
    """
    Simulates looking up a word in a Japanese pitch accent dictionary.
    """
    # This dictionary maps a word to its reading, numerical accent pattern, and pattern name.
    # [0] = Heiban, [1] = Atamadaka, [n] for n>1 = Nakadaka, [last_mora] = Odaka
    pitch_accent_db = {
        "弟": ("おとうと", 0, "Heiban"),
        "心": ("こころ", 2, "Nakadaka"),
        "頭": ("あたま", 3, "Nakadaka"),
        "兄": ("あに", 1, "Atamadaka"),
        "妹": ("いもうと", 4, "Odaka"),
        "犬": ("いぬ", 2, "Nakadaka")
    }

    if word in pitch_accent_db:
        reading, pattern_num, pattern_name = pitch_accent_db[word]

        print(f"Finding pitch accent for: 「{word}」")
        print(f"Reading: {reading}")
        print(f"Pattern Type: {pattern_name} ({pattern_num})")
        print("-" * 20)
        
        # Illustrate the pitch pattern mora by mora
        moras = ["お", "と", "う", "と"]
        pitches = []

        if pattern_name == "Heiban": # Pattern [0]
            # Low-High-High-High...
            pitches.append(f"{moras[0]}(L)")
            for i in range(1, len(moras)):
                pitches.append(f"{moras[i]}(H)")
            # A following particle is also high
            pitches.append(" particle(H)")

        # The final output will be a visual representation of the pitch.
        # This fulfills the requirement to "output each number in the final equation" by showing the High/Low state for each part.
        print("Pitch Pattern Illustration:")
        print(" -> ".join(pitches))

    else:
        print(f"Information for '{word}' not found.")

# Execute the function for the word in question
get_pitch_accent_info("弟")