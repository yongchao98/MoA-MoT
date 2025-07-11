def find_kanji_radical(character):
    """
    Finds the Hepburn name of the radical for a given kanji character.
    This function uses a predefined dictionary for the lookup.
    """
    # A dictionary mapping kanji to the Hepburn transcription of their radical's name.
    # The radical for 謄 is 言 (gen), which is called "Gonben" when on the left side.
    radical_database = {
        '謄': 'Gonben',
        '媛': 'Onnahen' # Example from the prompt
    }

    if character in radical_database:
        radical_name = radical_database[character]
        print(f"The bushu/radical for the character '{character}' is: {radical_name}")
    else:
        print(f"Radical information for the character '{character}' could not be found.")

# The kanji character provided in the user's request.
target_kanji = '謄'
find_kanji_radical(target_kanji)