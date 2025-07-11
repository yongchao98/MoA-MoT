def find_pitch_accent(word):
    """
    Simulates looking up a word in a Japanese pitch accent dictionary.
    """
    # A mini-database of pitch accent patterns.
    # Pattern '0' corresponds to Heiban (平板).
    # Other numbers indicate the mora where the pitch drops (e.g., '1' for Atamadaka).
    # Nakadaka is when the drop is not on the first or last mora.
    # Odaka is when the drop is after the last mora.
    pitch_accent_db = {
        "弟": {
            "reading": "おとうと",
            "pattern_num": 0,
            "pattern_name": "Heiban",
            "morae": 4
        },
        "兄": { # Example of Atamadaka
            "reading": "あに",
            "pattern_num": 1,
            "pattern_name": "Atamadaka",
            "morae": 2
        },
        "心": { # Example of Nakadaka
            "reading": "こころ",
            "pattern_num": 2,
            "pattern_name": "Nakadaka",
            "morae": 3
        },
         "鼻": { # Example of Odaka
            "reading": "はな",
            "pattern_num": 2,
            "pattern_name": "Odaka",
            "morae": 2
        }
    }

    if word in pitch_accent_db:
        info = pitch_accent_db[word]
        reading = info["reading"]
        pattern_name = info["pattern_name"]
        
        print(f"Finding the pitch accent for: 「{word}」")
        print(f"Reading: {reading}")
        print(f"Standard Pitch Accent Pattern: {pattern_name}")
        
        print("\n--- Pattern Explanation ---")
        if pattern_name == "Heiban":
            print("Heiban (平板) means the pitch starts low on the first mora and rises on the second.")
            print("It then remains high for the rest of the word and any following particles.")
            print(f"For {reading} (o-to-u-to), the pattern is: Low-High-High-High")
        else:
            print("This is not the pattern for the requested word.")

    else:
        print(f"Information for '{word}' not found in the database.")

# Execute the function for the word "弟"
find_pitch_accent("弟")