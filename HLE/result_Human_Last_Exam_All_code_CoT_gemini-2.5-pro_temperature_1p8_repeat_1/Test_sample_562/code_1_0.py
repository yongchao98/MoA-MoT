def find_radical_name(character):
    """
    Looks up the radical for a given Kanji character from a predefined dictionary.
    """
    kanji_database = {
        '謄': {
            'radical_char': '月',
            'hepburn_name': 'Nikuzuki',
            'description': 'When the radical 月 is on the left side, it is known as にくづき (nikuzuki).'
        },
        '媛': {
            'radical_char': '女',
            'hepburn_name': 'Onnahen',
            'description': 'The radical 女 (woman) on the left side.'
        }
    }

    if character in kanji_database:
        info = kanji_database[character]
        radical_char = info['radical_char']
        hepburn_name = info['hepburn_name']
        print(f"The bushu (radical) in the character '{character}' is '{radical_char}'.")
        print(f"The Hepburn transcription of its name is: {hepburn_name}")
    else:
        print(f"Information for the character '{character}' was not found.")

# The character we are analyzing
target_kanji = '謄'
find_radical_name(target_kanji)