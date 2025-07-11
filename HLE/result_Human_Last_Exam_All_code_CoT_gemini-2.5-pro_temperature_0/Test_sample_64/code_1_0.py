def find_pitch_accent(word):
    """
    Simulates looking up a word's pitch accent in a Japanese dictionary.
    """
    # A small database of pitch accent patterns.
    # L = Low pitch, H = High pitch
    # The notation [n] indicates the position of the pitch drop.
    # [0] = Heiban (平板): LHHHH... (no drop within the word)
    # [1] = Atamadaka (頭高): HLLLL... (drop after the 1st mora)
    # [2], [3], etc. = Nakadaka (中高): LHLL..., LHHLL... (drop in the middle)
    # Odaka (尾高) is when the drop occurs after the last mora, on the particle.
    pitch_accent_db = {
        "弟": {
            "reading": "おとうと (o-to-u-to)",
            "pattern_type": "Heiban",
            "pattern_number": 0,
            "explanation": "The pitch starts low, rises on the second mora, and stays high. (LHHH)"
        },
        "妹": {
            "reading": "いもうと (i-mo-u-to)",
            "pattern_type": "Heiban",
            "pattern_number": 0,
            "explanation": "The pitch starts low, rises on the second mora, and stays high. (LHHH)"
        },
        "頭": {
            "reading": "あたま (a-ta-ma)",
            "pattern_type": "Nakadaka",
            "pattern_number": 2,
            "explanation": "The pitch is Low-High-Low. (LHL)"
        },
        "命": {
            "reading": "いのち (i-no-chi)",
            "pattern_type": "Atamadaka",
            "pattern_number": 1,
            "explanation": "The pitch is High-Low-Low. (HLL)"
        },
        "花": {
            "reading": "はな (ha-na)",
            "pattern_type": "Odaka",
            "pattern_number": 2, # For a 2-mora word, a drop at [2] is Odaka
            "explanation": "The pitch is Low-High, and drops on a following particle. (LH-L)"
        }
    }

    if word in pitch_accent_db:
        info = pitch_accent_db[word]
        print(f"Searching for the pitch accent of 「{word}」...")
        print(f"Reading: {info['reading']}")
        print(f"Pitch Accent Type: {info['pattern_type']}")
        print(f"Explanation: {info['explanation']}")
        print("\nBased on the lookup, the correct pattern is Heiban.")
        print("This corresponds to Answer Choice A.")
    else:
        print(f"Word '{word}' not found in the database.")

# Execute the function for the target word
find_pitch_accent("弟")