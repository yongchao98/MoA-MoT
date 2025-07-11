def find_radical_name():
    """
    This function demonstrates finding the radical name for a given kanji.
    It uses a predefined dictionary to store the information.
    """
    # A dictionary mapping specific kanji to the Hepburn romanization of their radical's name.
    kanji_radical_map = {
        '謄': 'Gonben',
        '媛': 'Onnahen',
        '校': 'Kihen'
    }

    # The target character.
    character = '謄'

    # Retrieve and print the radical name.
    if character in kanji_radical_map:
        radical_name = kanji_radical_map[character]
        print(f"The bushu/radical for the character '{character}' is called: {radical_name}")
    else:
        print(f"Radical information for '{character}' not found.")

find_radical_name()