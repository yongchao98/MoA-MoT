def find_pitch_accent(word):
    """
    Simulates looking up the pitch accent of a Japanese word in a dictionary.
    """
    # Pitch accent database
    # Notation [0]: Heiban, [1]: Atamadaka, [n]: Nakadaka/Odaka
    pitch_accent_db = {
        '弟': {
            'reading': 'おとうと',
            'pattern_number': 0,
            'pattern_name': 'Heiban'
        },
        '箸': {
            'reading': 'はし',
            'pattern_number': 1,
            'pattern_name': 'Atamadaka'
        },
        '心': {
            'reading': 'こころ',
            'pattern_number': 2,
            'pattern_name': 'Nakadaka'
        },
        '鼻': {
            'reading': 'はな',
            'pattern_number': 2,
            'pattern_name': 'Odaka'
        }
    }

    if word in pitch_accent_db:
        info = pitch_accent_db[word]
        print(f"Searching for the pitch accent of: 「{word}」")
        print(f"Reading: {info['reading']}")
        print(f"The pitch accent pattern is [{info['pattern_number']}].")
        print(f"This corresponds to the '{info['pattern_name']}' pattern.")
    else:
        print(f"Information for 「{word}」 not found in this dictionary.")

# The word in question
target_word = '弟'
find_pitch_accent(target_word)