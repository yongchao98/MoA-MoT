def find_kanji_radical(character):
    """
    Finds the Hepburn transcription of the radical name for a given Kanji character.
    This function uses a predefined dictionary for lookup.
    """
    # A dictionary mapping Kanji characters to their radical's Hepburn name.
    # The name "Nikuzuki" is used for the '月' radical when it represents 'flesh' (肉).
    radical_database = {
        '謄': 'Nikuzuki',
        '媛': 'Onnahen', # Example from the prompt
        '海': 'Sanzui',
        '森': 'Ki'
    }

    if character in radical_database:
        radical_name = radical_database[character]
        print(f"The bushu/radical for the character '{character}' is '{radical_name}'.")
    else:
        print(f"Sorry, the radical for '{character}' is not in our database.")

# The specific character requested by the user.
target_character = '謄'
find_kanji_radical(target_character)