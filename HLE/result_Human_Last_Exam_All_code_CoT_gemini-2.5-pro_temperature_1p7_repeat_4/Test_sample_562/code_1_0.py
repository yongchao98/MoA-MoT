def find_kanji_radical_name(character):
    """
    Finds the Hepburn-romanized name of the radical for a given Kanji character.
    This function uses a predefined dictionary to simulate a radical database lookup.
    """
    kanji_radical_db = {
        '謄': 'Gonben', # Radical is 言 (gon) on the left (hen)
        '媛': 'Onnahen', # Radical is 女 (onna) on the left (hen)
        '語': 'Gonben', # Radical is 言 (gon) on the left (hen)
        '海': 'Sanzui', # Radical is 氵 (sanzui)
        '森': 'Ki'       # Radical is 木 (ki)
    }

    if character in kanji_radical_db:
        radical_name = kanji_radical_db[character]
        print(f"The radical for the character '{character}' is called '{radical_name}' in Hepburn transcription.")
    else:
        print(f"Sorry, the radical for '{character}' is not in our database.")

# The character from the user's request
target_kanji = '謄'

find_kanji_radical_name(target_kanji)