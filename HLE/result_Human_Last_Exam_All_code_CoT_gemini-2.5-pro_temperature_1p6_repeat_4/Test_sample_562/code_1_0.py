def find_radical_name(character):
    """
    Finds the Hepburn transcription of the radical name for a given Kanji character.
    This function uses a predefined dictionary for the lookup.
    """
    kanji_radical_database = {
        '謄': {
            'radical_char': '言',
            'radical_name_japanese': 'ごんべん',
            'radical_name_hepburn': 'Gonben'
        },
        '媛': {
            'radical_char': '女',
            'radical_name_japanese': 'おんなへん',
            'radical_name_hepburn': 'Onnahen'
        }
    }

    if character in kanji_radical_database:
        info = kanji_radical_database[character]
        print(f"The radical (bushu) in the character '{character}' is '{info['radical_char']}'.")
        print(f"Its common name is {info['radical_name_japanese']}.")
        print(f"The Hepburn transcription is: {info['radical_name_hepburn']}")
    else:
        print(f"Sorry, the radical for '{character}' is not in the database.")

# The character in question is '謄'.
find_radical_name('謄')
