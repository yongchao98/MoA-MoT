def get_pitch_accent(word):
    """
    A simple dictionary to store and retrieve pitch accent information.
    Accent patterns are represented numerically:
    [0] -> Heiban (平板): Low-High-High...
    [1] -> Atamadaka (頭高): High-Low-Low...
    [n] -> Nakadaka (中高) or Odaka (尾高): Accent on the nth mora.
    """
    accent_database = {
        "弟": {
            "reading": "おとうと",
            "pattern_name": "Heiban",
            "pattern_number": 0,
            "morae_count": 4,
            "explanation": "The pitch starts low, rises on the second mora, and remains high."
        }
    }

    if word in accent_database:
        info = accent_database[word]
        print(f"Word: 「{word}」")
        print(f"Reading: {info['reading']}")
        print(f"Pattern Name: {info['pattern_name']}")
        print("-" * 20)
        print("Pitch accent can be represented numerically by the position of the accent kernel.")
        print("For Heiban, the accent kernel is on position 0, meaning there is no downstep in pitch within the word.")
        print(f"Numerical Pattern: [{info['pattern_number']}]")
        print("\nPronunciation Pattern (L = Low, H = High):")
        print("Mora 1 (お): L")
        print("Mora 2 (と): H")
        print("Mora 3 (う): H")
        print("Mora 4 (と): H")
    else:
        print(f"Pitch accent for '{word}' not found.")

# Find the pitch accent for 「弟」
get_pitch_accent("弟")